#!/usr/bin/env python

##   Copyright 2014 Bioinformatics Core, University of Michigan
##
##   Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##
##       http://www.apache.org/licenses/LICENSE-2.0
##
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

#pylint:disable=line-too-long, no-member

from __future__ import absolute_import

import argparse
from collections import OrderedDict
from datetime import datetime
import errno
import logging
import os
import re
import shlex
import sys

from colour import Color
import numpy as np
import pandas as pd

from generollup import __version__

_LOG_DATE_FORMAT = '%Y-%m-%d %H:%M:%S'
_LOG_LINE_FORMAT = '%(asctime)s|%(levelname)s|%(message)s'

# input columns
#_GENE_SYMBOL = 'Gene Names'
#_SAMPLE_GENOTYPE_COLUMN_REGEX = r'(.*) 0/1 Genotypes \(GT\)'
#_DBNSFP_COLUMN = 'N of 6 Predicted Damaging'
#_EFFECT_COLUMN = 'Effect (Combined)'

# output columns
_GENE_SYMBOL_OUTPUT_NAME = 'gene_symbol'
_LOCI_COUNT = 'distinct_variant_loci'
_SAMPLE_COUNT = 'total_affected_samples'
_MUTATION_COUNT = 'total_mutations'


_DESIRED_ANNOTATIONS = set()


def _is_variant(genotype):
    return sum(int(x) for x in genotype if x.isdigit()) > 0

def _is_not_variant(genotype):
    return not _is_variant(genotype)

class RollupException(Exception):
    """Base class for all run-time exceptions in this module."""
    def __init__(self, msg, *args):
        #pylint: disable=star-args
        error_msg = msg.format(*[str(i) for i in args])
        super(RollupException, self).__init__(error_msg)

class UsageError(RollupException):
    """Raised for malformed command or invalid arguments."""
    def __init__(self, msg, *args):
        super(UsageError, self).__init__(msg, *args)

def _count_values(df, column, cutoff=20):
    counts = df[column].value_counts().to_dict()
    counts = sorted(counts.items(), key=lambda x: -1 * x[1])
    num_counts = len(counts)
    msgs = [(column + ' value', str(counts[i][0]), str(counts[i][1])) for i in range(0, min(cutoff, num_counts))]
    omitted = ''
    if num_counts > cutoff:
        omitted = '|values omitted from log|{}'.format(num_counts - cutoff)
    msgs.append((column, 'distinct values', '{}{}'.format(num_counts, omitted),))
    return msgs

class dbNSFP(object):
    #pylint: disable=invalid-name, too-few-public-methods
    def __init__(self, format_rules, args):
        self.name = 'dbNSFP_annotation'
        self.args = args
        self.column_label = 'damaging_votes'
        self.damaging_column = _DBNSFP_COLUMN
        self.damaging_rank_column = 'overall_damaging_rank'
        self.damaging_total = 'damaging_total'
        self.format_rules = format_rules

    def validation(self, df):
        msgs = [('damaging_column', self.damaging_column)]
        msgs.extend(_count_values(df, self.damaging_column))
        return msgs

    @staticmethod
    def _to_rounded_string(cell_value):
        try:
            return re.findall(r'-?\d+', str(cell_value))[0]
        except IndexError:
            return ''

    def summarize(self, initial_df):
        condensed_df = self._remove_unnecessary_columns(initial_df)
        validated_df = self._remove_invalid_rows(condensed_df)

        damaging_votes_df = self._build_gene_sample_damaging_votes(validated_df)
        damaging_ranks_df = self._build_damaging_ranks(damaging_votes_df)

        frames = [damaging_ranks_df, damaging_votes_df]
        data_df = pd.concat(frames, axis=1)

        data_df = data_df.applymap(self._to_rounded_string)
        data_df = data_df.replace(to_replace=0,value='')

        style_dfs = []
        for format_rule in self.format_rules:
            style_dfs.append(format_rule.format(data_df))

        return data_df, style_dfs

    def _remove_unnecessary_columns(self, initial_df):
        sample_cols = initial_df.filter(regex=self.args.sample_genotype_column_regex)
        required_columns = list(sample_cols.columns)
        required_columns.extend([self.args.gene_column_name,
                                 self.damaging_column])

        condensed_df = pd.DataFrame()
        for col in initial_df.columns:
            if col in required_columns:
                condensed_df[col] = initial_df[col]

        return condensed_df

    def _remove_invalid_rows(self, initial_df):
        initial_df.fillna('.', inplace=True)
        initial_df = initial_df[initial_df[self.damaging_column] != '.']

        #temporarily reset warning thresholds to prevent spurious warning in pandas
        default = pd.options.mode.chained_assignment
        try:
            pd.options.mode.chained_assignment = None
            initial_df[self.damaging_column] = initial_df.loc[:,self.damaging_column].apply(lambda x: int('0' + self._to_rounded_string(x)))
            initial_df = initial_df[initial_df[self.damaging_column] > 0]
            initial_df[self.damaging_column] = initial_df[self.damaging_column].apply(str)
        finally:
            pd.options.mode.chained_assignment = default

        return initial_df

    def _build_gene_sample_damaging_votes(self, initial_df):
        sample_cols = initial_df.filter(regex=self.args.sample_genotype_column_regex).columns.values
        initial_df = initial_df.applymap(str)

        for sample in sample_cols:
            initial_df[sample][initial_df[sample].apply(_is_not_variant)] = '.'

            initial_df.fillna('.', inplace=True)
            initial_df[sample][initial_df[sample] != '.'] = initial_df[self.damaging_column]

            samp_suffix = re.match(self.args.sample_genotype_column_regex, sample)[1]
            sample_name = '|'.join([self.name,self.column_label,samp_suffix])
            initial_df[sample] = initial_df[sample].apply(lambda x: x if x != '.' else np.nan)
            initial_df[sample_name] = initial_df[sample].apply(float)
            del initial_df[sample]

        damaging_votes_df = initial_df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
        damaging_votes_df = damaging_votes_df.applymap(lambda x: x if x != '.' else 0)
        damaging_votes_df = damaging_votes_df.groupby(self.args.gene_column_name).sum()
        damaging_votes_df = damaging_votes_df.applymap(lambda x: x if x != 0 else '.')
        del damaging_votes_df[self.damaging_column]

        return damaging_votes_df

    def _build_damaging_ranks(self,damaging_votes_df):
        damaging_rank_name = '|'.join([self.name, self.damaging_rank_column])
        total_damaging_name = '|'.join([self.name, self.damaging_total])

        damaging_ranks_df = pd.DataFrame()
        damaging_ranks_df[total_damaging_name] = damaging_votes_df.apply(pd.to_numeric, errors='coerce').sum(axis=1, skipna=True)
        damaging_ranks_df[total_damaging_name].fillna(0, inplace=True)

        damaging_ranks_df[damaging_rank_name] = damaging_ranks_df[total_damaging_name].rank(ascending=False,
                                                                                            method='min')
        damaging_ranks_df[damaging_rank_name] = damaging_ranks_df[damaging_rank_name].apply(np.floor)
        damaging_ranks_df = damaging_ranks_df.loc[:,[damaging_rank_name, total_damaging_name]]

        return damaging_ranks_df

class Effect(object):
    #pylint: disable=invalid-name, too-few-public-methods
    # _RANK_SCORES = {'LoF': 10**9,
    #                 'Missense': 10**6,
    #                 'Other': 10**3,
    #                 'Missing': 10**0}
    # _RANK_ABBREVS = {'LoF': 'h',
    #                  'Missense': 'm',
    #                  'Other': 'l',
    #                  'Missing': 'x',
    #                  '.': ''}

    def __init__(self, format_rules, args):
        self.name = 'effect_annotation'
        self.args = args
        self.column_label = 'effect'
        self.effect_column = _EFFECT_COLUMN
        self.effect_rank_column = 'overall_effect_rank'
        self.effect_score_column = 'overall_effect_score'
        self.effect_category = 'effect_category'
        self.format_rules = format_rules
        #TODO: cgates: "If guard" is ugly hack to allow the no-arg constructors in styling code
        if args:
            self._rank_scores = {args.effect_column_values[0]: 10**9,
                                 args.effect_column_values[1]: 10**6,
                                 args.effect_column_values[2]: 10**3,
                                 args.effect_column_values[3]: 10**0}
            self._rank_abbrevs = {args.effect_column_values[0]: 'h',
                                  args.effect_column_values[1]: 'm',
                                  args.effect_column_values[2]: 'l',
                                  args.effect_column_values[3]: 'x',
                                  '.': ''}

    def validation(self, df):
        msgs = [('effect_column', self.effect_column)]
        msgs.extend(_count_values(df, self.effect_column))
        return msgs

    def summarize(self, initial_df):
        condensed_df = self._remove_unnecessary_columns(initial_df)
        validated_df = self._remove_invalid_rows(condensed_df)

        scored_df = self._calculate_score(validated_df)
        category_df = self._get_effect_category_counts(scored_df)
        ranked_df = self._calculate_rank(category_df)
        data_df = self._change_col_order(ranked_df)

        style_dfs = []
        for format_rule in self.format_rules:
            style_dfs.append(format_rule.format(data_df))

        return data_df, style_dfs

    def _remove_unnecessary_columns(self, initial_df):
        sample_cols = initial_df.filter(regex=self.args.sample_genotype_column_regex)
        required_columns = list(sample_cols.columns)
        required_columns.extend([self.args.gene_column_name,
                                 self.effect_column])

        condensed_df = pd.DataFrame()
        for col in initial_df.columns:
            if col in required_columns:
                condensed_df[col] = initial_df[col]

        return condensed_df

    def _remove_invalid_rows(self, initial_df):
        initial_df.fillna('.', inplace=True)
        initial_df = initial_df[initial_df[self.effect_column] != '.']

        return initial_df

    def _calculate_score(self, initial_df):
        sample_regex = self.args.sample_genotype_column_regex
        sample_cols = initial_df.filter(regex=sample_regex).columns.values
        initial_df = initial_df.applymap(str)

        scored_df = pd.DataFrame()
        scored_df[_GENE_SYMBOL_OUTPUT_NAME] = initial_df[self.args.gene_column_name]

        for sample in sample_cols:
            #set sample columns equal to effect column value
            initial_df[sample][initial_df[sample].apply(_is_not_variant)] = '.'

            initial_df.fillna('.', inplace=True)
            initial_df[sample][initial_df[sample] != '.'] = initial_df[self.effect_column]

            scored_df[sample + '_SCORE'] = initial_df[sample].map(self._rank_scores)
            initial_df[sample] = initial_df[sample].map(self._rank_abbrevs)

        initial_df.fillna('', inplace=True)
        scored_df.fillna(0, inplace=True)

        score = scored_df.groupby(_GENE_SYMBOL_OUTPUT_NAME).sum().apply(sum,1)
        score = score.apply(int)
        grouped_df = initial_df.groupby(self.args.gene_column_name).sum()
        grouped_df = grouped_df[grouped_df.index != '.']
        grouped_df = grouped_df.applymap(lambda x: ''.join(sorted(x)))

        if self.effect_column in grouped_df.columns.values:
            del grouped_df[self.effect_column]

        effect_score_name = '|'.join([self.name, self.effect_score_column])
        grouped_df[effect_score_name] = score
        grouped_df = grouped_df.applymap(str)

        return grouped_df

    def _get_effect_category_counts(self, initial_df):
        def _count_category(abbrev):
            sample_regex = self.args.sample_genotype_column_regex
            sample_cols = initial_df.filter(regex=sample_regex).columns.values
            temp_df = pd.DataFrame()
            for sample_col in sample_cols:
                temp_df[sample_col] = initial_df[sample_col].apply(lambda x: len([i for i in x if i == abbrev]))
            return temp_df.apply(sum, 1)

        effect_category_name = '|'.join([self.name, self.effect_category])
        category_df = initial_df.copy()
        category_df[effect_category_name + '|LoF'] = _count_category('h')
        category_df[effect_category_name + '|Missense'] = _count_category('m')
        category_df[effect_category_name + '|Other'] = _count_category('l')
        category_df[effect_category_name + '|Missing'] = _count_category('x')

        return category_df

    def _rename_sample_columns(self, initial_df):
        sample_regex = self.args.sample_genotype_column_regex
        sample_cols = initial_df.filter(regex=sample_regex).columns.values
        for sample in sample_cols:
            sample_suffix = re.match(sample_regex, sample)[1]
            sample_name = '|'.join([self.name,
                                    self.column_label,
                                    sample_suffix])
            initial_df[sample_name] = initial_df[sample]
            del initial_df[sample]

        return initial_df

    def _calculate_rank(self, initial_df):
        effect_rank_name = '|'.join([self.name, self.effect_rank_column])
        effect_score_name = '|'.join([self.name, self.effect_score_column])

        category_df = initial_df
        category_df[effect_score_name] = category_df[effect_score_name].apply(int)
        category_df[effect_rank_name] = category_df[effect_score_name].rank(ascending=False,
                                                                            method='min')
        category_df[effect_rank_name] = category_df[effect_rank_name].apply(int)
        category_df = category_df.applymap(str)

        return self._rename_sample_columns(category_df)

    def _change_col_order(self, initial_df):
        effect_rank_name = '|'.join([self.name, self.effect_rank_column])
        effect_score_name = '|'.join([self.name, self.effect_score_column])

        rank = []
        score = []
        category = []
        effect = []

        for col in initial_df.columns.values:
            if col == effect_rank_name:
                rank.append(col)
            elif col == effect_score_name:
                score.append(col)
            else:
                regex = r'^{}\|{}\|.*'.format(self.name, self.column_label)
                if re.match(regex, col):
                    effect.append(col)
                else:
                    category.append(col)

        ordered_names = rank + score + category + effect
        initial_df = initial_df[ordered_names]

        return initial_df


class SummaryColumns(object):
    def __init__(self, args):
        self.name = 'summary_annotation'
        self.args = args

    def validation(self, df):
        msgs = [('sample_genotype_column_regex', self.args.sample_genotype_column_regex),
                ('gene_column_name', self.args.gene_column_name),
               ]
        sample_df = df.filter(regex=self.args.sample_genotype_column_regex,
                              axis=1)
        unique_gt = sorted(np.unique(sample_df.values))
        variant_gt = [('variant genotypes',
                       str(list(filter(_is_variant, unique_gt)))),
                      ('non-varaint genotypes',
                       str(list(filter(lambda x: not _is_variant(x), unique_gt)))),
                     ]
        msgs.extend(variant_gt)
        return msgs

    def summarize(self, initial_df):
        sample_df = initial_df.filter(regex=self.args.sample_genotype_column_regex).copy()
        sample_df[_GENE_SYMBOL_OUTPUT_NAME] = initial_df[self.args.gene_column_name]

        data_df = pd.DataFrame()
        data_df[_SAMPLE_COUNT] = self.calculate_affected_samples(sample_df)
        data_df[_LOCI_COUNT] = self.calculate_variant_loci(sample_df)
        data_df[_MUTATION_COUNT] = self.calculate_total_mutations(sample_df)
        data_df = data_df.applymap(str)

        style_dfs = []
        style_df = data_df.copy()
        style_dfs.append(style_df.applymap(lambda x: ''))

        return data_df, style_dfs

    def calculate_affected_samples(self, sample_df):
        sample_df = sample_df.set_index(_GENE_SYMBOL_OUTPUT_NAME).applymap(_is_variant)
        return sample_df.groupby(sample_df.index).any().sum(axis=1)

    def calculate_variant_loci(self, sample_df):
        sample_df = sample_df.set_index('gene_symbol').applymap(_is_variant)
        return sample_df.sum(axis=1).apply(lambda x: 1 if x>0 else 0).groupby(level=0).sum()

    def calculate_total_mutations(self, sample_df):
        sample_df = sample_df.set_index(_GENE_SYMBOL_OUTPUT_NAME).applymap(_is_variant)
        return sample_df.groupby(_GENE_SYMBOL_OUTPUT_NAME).sum().sum(axis=1)


class EffectFormatRule(object):
    #pylint: disable=too-few-public-methods
    _COLOR_MAP = {'h': '#003366',
                  'm': '#3377B9',
                  'l': '#91C4E8',
                  'x': '#CCCCCC'}
    _FONT_COLOR = '#000000'

    def __init__(self):
        self.rank_column = 'effect|overall_effect_rank'
        self.range_rule = list(Color('white').range_to(Color('blue'), 101))

    def format(self, data_df):
        def _determine_format(cell_value):
            for letter in EffectFormatRule._COLOR_MAP.keys():
                if letter in cell_value:
                    return letter
            return ''

        format_df = data_df.applymap(str)
        format_df = format_df.applymap(_determine_format)
        format_df = format_df.applymap(self._style)

        return format_df

    @staticmethod
    def _style(cell_value):
        if len(cell_value) < 1:
            return np.nan

        else:
            styled_cell = {'font_size': '4',
                           'bg_color': EffectFormatRule\
                                               ._COLOR_MAP[cell_value],
                           'font_color': EffectFormatRule\
                                         ._COLOR_MAP[cell_value]}
            return styled_cell


class dbNSFPFormatRule(object):
    #pylint: disable=invalid-name,too-few-public-methods

    _FONT_COLOR = '#000000'

    def __init__(self):
        self.rank_column = 'dbNSFP_annotation|overall_damaging_rank'
        self.total_column = 'dbNSFP_annotation|damaging_total'
        self.range_rule = list(Color('white').range_to(Color('orange'), 101))

    def format(self, data_df):
        def _determine_format(cell_value):
            normalized = 100*(cell_value-min_value)/(max_value-min_value)
            if np.isnan(normalized):
                return cell_value
            else:
                return int(normalized)

        format_df = data_df.drop(self.rank_column, 1)
        format_df = format_df.drop(self.total_column, 1)

        format_df = format_df.replace('', np.nan)
        format_df = format_df.apply(pd.to_numeric)

        min_value = min(format_df.min(skipna=True))
        max_value = max(format_df.max(skipna=True))

        format_df = format_df.applymap(_determine_format)
        format_df[self.rank_column] = pd.Series()
        format_df[self.total_column] = pd.Series()

        format_df = format_df.applymap(self._style)

        return format_df

    def _style(self, cell_value):
        if np.isnan(cell_value):
            return np.nan
        else:
            cell_value = int(cell_value)
            color = self.range_rule[cell_value]
            styled_cell = {'font_size': '12',
                           'bg_color': str(color),
                           'font_color': dbNSFPFormatRule._FONT_COLOR}
            return styled_cell


class RankFormatRule(object):
    #pylint: disable=invalid-name, too-few-public-methods

    _FONT_COLOR = '#000000'

    def __init__(self):
        self.range_rule = list(Color('red').range_to(Color('white'), 101))

    def format(self, data_df):
        def _determine_format(cell_value):
            if len(str(cell_value)) > 0:
                cell_value = int(cell_value)
                normalized = 100*(cell_value-min_value)/(max_value-min_value)
                return str(int(normalized))
            else:
                return ''

        format_df = data_df.copy()

        rank_columns = {}

        for col in format_df.columns.values:
            if re.search(dbNSFP('','').damaging_rank_column, col) or \
                re.search(Effect('','').effect_rank_column, col):
                min_value = int(format_df[col].min(skipna=True))
                max_value = int(format_df[col].max(skipna=True))

                if max_value - min_value > 0:
                    format_df.fillna('', inplace=True)
                    format_df[col] = format_df[col].apply(_determine_format)

                    rank_columns[col] = format_df[col]

                else:
                    format_df[col] = format_df[col].apply(lambda x: '')

            else:
                format_df[col] = format_df[col].apply(lambda x: '')

        try:
            max_range_value = int(max(format_df.max(1)))
            self._set_max_range_rule(max_range_value)
        except ValueError:
            pass

        format_df = pd.DataFrame(data=None,
                                 columns=format_df.columns,
                                 index=format_df.index)

        #pylint: disable=unnecessary-lambda
        for rank_col_name, rank_col in rank_columns.items():
            format_df[rank_col_name] = rank_col.apply(lambda x: self._style(x))

        return format_df

    def _set_max_range_rule(self, value):
        value = value + 1
        self.range_rule = list(Color('#E27171').range_to(Color('white'), value))

    def _style(self, cell_value):
        if len(str(cell_value)) > 0:
            cell_value = int(cell_value)
            color = self.range_rule[cell_value]
            styled_cell = {'font_size': '12',
                           'bg_color': str(color),
                           'font_color': RankFormatRule._FONT_COLOR}
            return styled_cell
        return np.nan


def _create_df(input_file, args):
    initial_df = pd.read_csv(input_file, sep='\t', header=0, dtype='str')
    _info('read {} ({} variant loci x {} columns)', input_file, *initial_df.shape)

    #TODO: cgates: as implemented, you must do this, but consider an refactor which would remove this side-effect
    if args.gene_column_index:
        args.gene_column_name = initial_df.columns.values[args.gene_column_index-1]
        global _GENE_SYMBOL
        _GENE_SYMBOL = args.gene_column_name
        _info('using input column index={}, name=[{}] for gene name',
              args.gene_column_index,
              args.gene_column_name)
    _validate_df(initial_df, args)
    initial_df = initial_df[pd.notnull(initial_df[_GENE_SYMBOL])]

    return initial_df

def _validate_df(initial_df, args):
    required_columns = {}

    required_columns['gene_column_name'] = args.gene_column_name

    if args.dbnsfp_column_name:
        _DESIRED_ANNOTATIONS.add(args.dbnsfp_column_name)
        required_columns['dbnsfp_column_name'] = args.dbnsfp_column_name

    if args.effect_column_name:
        _DESIRED_ANNOTATIONS.add(args.effect_column_name)
        required_columns['effect_column_name'] = args.effect_column_name

    missing_columns = {}
    header = initial_df.columns.values
    for option, column_name in required_columns.items():
        if column_name not in header:
            missing_columns[option] = column_name

    msg = ('Some options referenced column names not found in input file: {}. '
           'Review inputs and try again.')\
           .format(', '.join(['{}:{}'.format(k,v) for k,v in missing_columns.items()]))
    if missing_columns:
        raise UsageError(msg)


    sample_df = initial_df.filter(regex=args.sample_genotype_column_regex, axis=1)
    if not len(sample_df.columns):
        msg = ('Cannot determine sample genotype columns with supplied regex ({}). '
               'Review input and try again.').format(args.sample_genotype_column_regex)
        raise UsageError(msg)

    _info('{} variant loci x {} samples', *sample_df.shape)

def _combine_style_dfs(dfs):
    df1, df2 = dfs
    original_columns = set(df1.columns.values).union(set(df2.columns.values))

    cleaned_df1 = df1.dropna(1, how='all')
    cleaned_df2 = df2.dropna(1, how='all')

    combined_df = cleaned_df1.join(cleaned_df2, how='outer')
    for column in original_columns:
        if column not in list(combined_df.columns.values):
            combined_df[column] = pd.Series()

    combined_df = combined_df.fillna('')

    return combined_df

def _combine_dfs(dfs):
    summary_df = list(dfs.values())[0]
    annotation_dfs = list(dfs.values())[1:]

    summary_genes = set(summary_df.index.values)

    valid_genes = set()
    for annotation_df in annotation_dfs:
        df_genes = set(annotation_df.index.values)
        valid_genes = valid_genes.union(df_genes)

    invalid_genes = summary_genes.difference(valid_genes)
    summary_df = summary_df.drop(list(invalid_genes))

    combined_df = summary_df
    for annotation_df in annotation_dfs:
        combined_df = combined_df.join(annotation_df, how='outer')

    combined_df = combined_df.fillna('')
    combined_df.index.names = [_GENE_SYMBOL_OUTPUT_NAME]

    combined_df = combined_df[combined_df.index != '.']

    return combined_df

def _header_formats():
    #pylint: disable=duplicate-key
    dbnsfp = dbNSFP('','')
    effect = Effect('','')
    header_formats = {_GENE_SYMBOL_OUTPUT_NAME: {'align': 'center',
                                                 'bold': True,
                                                 'border': True,
                                                 'align': 'center'},
                      dbnsfp.damaging_rank_column: {'bg_color': '#F8A049',
                                                    'font_color': 'white',
                                                    'border': True,
                                                    'rotation': 90,
                                                    'align': 'center'},
                      dbnsfp.damaging_total: {'bg_color': '#FBC692',
                                              'border': True,
                                              'rotation': 90,
                                              'align': 'center'},
                      dbnsfp.column_label: {'bg_color': '#FBC692',
                                            'border': True,
                                            'rotation': 90,
                                            'align': 'center'},
                      effect.effect_rank_column: {'bg_color': '#4775A3',
                                                  'font_color': 'white',
                                                  'border': True,
                                                  'rotation': 90,
                                                  'align': 'center'},
                      effect.effect_score_column: {'bg_color': '#6C91B5',
                                                   'font_color': 'white',
                                                   'border': True,
                                                   'rotation': 90,
                                                   'align': 'center'},
                      effect.effect_category: {'bg_color': '#6C91B5',
                                               'font_color': 'white',
                                               'border': True,
                                               'rotation': 90,
                                               'align': 'center'},
                      effect.column_label + r'\|': {'bg_color': '#99C1D7',
                                                    'border': True,
                                                    'rotation': 90,
                                                    'align': 'center'},
                      _LOCI_COUNT: {'bg_color': '#F7F4FF',
                                    'border': True,
                                    'rotation': 90,
                                    'align': 'center'},
                      _SAMPLE_COUNT: {'bg_color': '#F7F4FF',
                                      'border': True,
                                      'rotation': 90,
                                      'align': 'center'},
                      _MUTATION_COUNT: {'bg_color': '#F7F4FF',
                                        'border': True,
                                        'rotation': 90,
                                        'align': 'center'}}
    return header_formats

def _translate_to_excel(data_df, style_df, writer):
    gene_symbols = data_df[_GENE_SYMBOL_OUTPUT_NAME].values
    data_df = data_df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
    data_df[_GENE_SYMBOL_OUTPUT_NAME] = gene_symbols

    worksheet_name = 'gene_rollup'
    data_df.to_excel(writer, sheet_name=worksheet_name, index=False)

    workbook = writer.book
    workbook.strings_to_numbers = True
    worksheet = writer.sheets[worksheet_name]
    worksheet.strings_to_numbers = True

    worksheet.set_column(0, 0, 20)
    worksheet.set_column(1, len(list(data_df.columns.values)), 5)

    for i, column in enumerate(data_df.columns.values):
        for key in _header_formats():
            if re.search(dbNSFP('','').damaging_total, column) or \
                re.search(Effect('','').effect_score_column, column):
                worksheet.set_column(i, i, None, None, {'hidden': 1})

            if re.search(key, column):
                cell_format = workbook.add_format(_header_formats()[key])
                worksheet.write(0, i, column, cell_format)
                break

    _info('translating to Excel format')
    for i, (row, dummy) in enumerate(data_df.iterrows()):
        for j, (column, dummy) in enumerate(data_df.iteritems()):
            style = style_df.loc[row, column]
            if style:
                cell_format = workbook.add_format(style)
                worksheet.write(i + 1, j, data_df.loc[row, column], cell_format)


    writer.save()

def _sort_by_dbnsfp_rank(initial_df):
    dbnsfp = dbNSFP('','')
    effect = Effect('','')
    dbnsfp_rank_col = '|'.join([dbnsfp.name, dbnsfp.damaging_rank_column])
    effect_rank_col = '|'.join([effect.name, effect.effect_rank_column])

    sort_order = []
    try:
        initial_df[dbnsfp_rank_col] = initial_df[dbnsfp_rank_col].replace('',
                                                                          np.nan)
        initial_df[dbnsfp_rank_col] = initial_df[dbnsfp_rank_col].apply(float)
        sort_order.append(dbnsfp_rank_col)
    except KeyError:
        pass

    try:
        initial_df[effect_rank_col] = initial_df[effect_rank_col].replace('',
                                                                          np.nan)
        initial_df[effect_rank_col] = initial_df[effect_rank_col].apply(float)
        sort_order.append(effect_rank_col)
    except KeyError:
        pass

    sort_order.append(_GENE_SYMBOL_OUTPUT_NAME)
    sorted_df = initial_df.sort_values(by=sort_order)

    sorted_df.fillna('', inplace=True)

    return sorted_df

def _reset_style_gene_values(combined_style_df):
    original_values = combined_style_df[_GENE_SYMBOL_OUTPUT_NAME]
    new_values = [{} for dummy in original_values]
    combined_style_df[_GENE_SYMBOL_OUTPUT_NAME] = new_values

    return combined_style_df

def _rollup(args):
    _info('gene_rollup ({}) begins', __version__)
    _info('command | {}', ' '.join(map(shlex.quote, sys.argv[1:])))
    initial_df = _create_df(args.input_file, args)

    annotations = [SummaryColumns(args)]
    if _DBNSFP_COLUMN in _DESIRED_ANNOTATIONS:
        annotations.append(dbNSFP([dbNSFPFormatRule(), RankFormatRule()], args))
    if _EFFECT_COLUMN in _DESIRED_ANNOTATIONS:
        annotations.append(Effect([EffectFormatRule(), RankFormatRule()], args))

    annotation_dfs = OrderedDict()
    all_style_dfs = OrderedDict()

    for annotation in annotations:
        for msgs in annotation.validation(initial_df):
            _info('annotation|{}|validation|{}', annotation.name, '|'.join(msgs))
        summarized_df, style_dfs = annotation.summarize(initial_df)
        annotation_dfs[annotation.name] = summarized_df
        if len(style_dfs) == 2:
            combined_style_df = _combine_style_dfs(style_dfs)
            all_style_dfs[annotation.name] = combined_style_df
        elif len(style_dfs) == 1:
            all_style_dfs[annotation.name] = style_dfs[0]

    combined_df = _combine_dfs(annotation_dfs)
    combined_df = combined_df.reset_index()
    sorted_df = _sort_by_dbnsfp_rank(combined_df)

    combined_style_df = _combine_dfs(all_style_dfs)
    combined_style_df = combined_style_df.reset_index()
    altered_style_df = _reset_style_gene_values(combined_style_df)

    _info('writing to output file [{}] ({} x {})',
          args.output_file,
          *sorted_df.shape)

    if args.tsv:
        sorted_df.to_csv(args.output_file, sep='\t', index=False)
    else:
        try:
            writer = pd.ExcelWriter(args.output_file, engine='xlsxwriter')
            _translate_to_excel(sorted_df, altered_style_df, writer)

        except ValueError:
            msg = ('Unable to write [{}] to an Excel file. Review inputs and '
                   'try again.').format(args.output_file)
            raise RollupException(msg)

    _info('done')

def _parse_args(args):
    parser = argparse.ArgumentParser(\
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=\
'''Aggregates a table of annotated variants x samples into a table of genes x
samples. The rollup expects an input file where each variant locus is
represented once and each locus is assigned to a single gene; for example, your
input could be the highest impact alt allele in the highest impact transcript of
a gene. That said, rollup will aggregate a locus assigned to multiple genes or
transcripts; these cases may simply require more careful interpretation.

Based on command line args will optionally summarize several annotation types
for each impacted gene:
- summary effects (num affected samples, num distinct loci, num total variants)
- categorical effects (ala SnpEff (HIGH,MODERATE,LOW,MODIFIER) or
  VarSeq (LoF,Misense,Other,Missing))
- numeric impact (ala count of dbSNP damaging votes)

Note that the input file is not filtered by rows, so apply upstream filters
as necessary (e.g. VCF FILTER=PASS or filter to gene regions).

By default emits styled Excel workbook, but also supports TSV output.

See doc for detailed explanation on how various output columns are derived.''',
        epilog='See https://github.com/umich-brcf-bioinf/GeneRollup for more info.',
        )
    parser.add_argument('input_file',
                        help=('An annotated tab-separated file of variant '
                              'loci x samples.'))
    parser.add_argument('output_file',
                        help=('By default, an XLSX file of genes x samples. '
                              'Can be a TSV file with --tsv'))
    parser.add_argument('--sample_genotype_column_regex',
                        required=True,
                        #default=_SAMPLE_GENOTYPE_COLUMN_REGEX,
                        help=('Regex used to define the sample genotype columns in the '
                              'input file. Any non-zero number in this column will flag '
                              'the loci as a variant; so "0/1", "1|0", "2" are all variants '
                              'but "0/0", "./.", ".", "0" are not. Also, the first group in '
                              'the regex will be used as the sample name in the output columns.'))
    parser.add_argument('--tsv',
                        action='store_true',
                        default=False,
                        help='Write to a tsv file rather than an xlsx file')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--gene_column_name',
                        #default=_GENE_SYMBOL,
                        help=('Gene column name in input file. '
                              'Each field value of this column should contain '
                              'exactly one gene symbol. (Note that all variants '
                              'with gene="." will be grouped together, so it''s '
                              'best to limit results to exclude intergenic '
                              'regions.)')
                        )
    group.add_argument('--gene_column_index',
                        type=int,
                        help=('One-based index of gene symbol column in the '
                              'input file; useful if the input column names '
                              'are not unique.'))
    parser.add_argument('--dbnsfp_column_name',
                        #default=_DBNSFP_COLUMN,
                        help=('dbNSFP damaging count column in the input file. '
                              'Each field value of this column must contain '
                              'at least one number, the first number is '
                              'is interpreted at the number of damaging votes ('
                              'i.e. "3 of 6" -> 3 annotation sources consider '
                              'this variant damaging). Note that "3", '
                              '"3 of 6" or "0" are all valid, but "none", ".", '
                              'or "" would raise an error.')
        )
    parser.add_argument('--effect_column_name',
                        #default=_EFFECT_COLUMN,
                        help=('Name of effect column in the input file which '
                              'contains categorial impact values (LoF, '
                              'Missense, etc.)'))
    parser.add_argument('--effect_column_values',
                        help=('Comma-separated list of exactly four '
                              'distinct effect values. Values should be listed '
                              'from highest to lowest impact, e.g. '
                              '"Lof,Missense,Other,Missing".'))
    parser.add_argument('--log_file',
                        help='defaults to <output_file>.log')
    parser.add_argument('--version',
                        '-V',
                        action='version',
                        version=__version__)
    args = parser.parse_args(args)
    if not args.log_file:
        args.log_file = args.output_file + '.log'
    if args.effect_column_name and args.effect_column_values:
        args.effect_column_values = args.effect_column_values.split(',')
    elif args.effect_column_name or args.effect_column_values:
        parser.error('effect annotation requires both effect_column_name and effect_column_values')
    return args

def _format(message, args):
    try:
        log_message = message.format(*[str(i) for i in args])
    except IndexError as err:
        log_message = ('Malformed log message ({}: {})'
                       '|{}|{}').format(type(err).__name__,
                                        err,
                                        message,
                                        [str(i) for i in args])
    return log_message

def _info(message, *args):
    logging.info(_format(message, args))

def main():
    args = _parse_args(sys.argv[1:])

    logging.basicConfig(format=_LOG_LINE_FORMAT,
                        level='DEBUG',
                        datefmt=_LOG_DATE_FORMAT,
                        filename=args.log_file)
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    formatter = logging.Formatter(_LOG_LINE_FORMAT, _LOG_DATE_FORMAT)
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    #TODO: cgates: This use of global is a dumb pattern. Let's get rid of it.
    global _GENE_SYMBOL
    _GENE_SYMBOL = args.gene_column_name
    global _DBNSFP_COLUMN
    _DBNSFP_COLUMN = args.dbnsfp_column_name
    global _EFFECT_COLUMN
    _EFFECT_COLUMN = args.effect_column_name

    #TODO: cgates: How about we fix the chained assigments and get rid of this warning suppression.
    pd.set_option('mode.chained_assignment', None)
    try:
        _rollup(args)
    except UsageError as e:
        print('Usage Error: ' + str(e))
if __name__ == '__main__':
    main()
