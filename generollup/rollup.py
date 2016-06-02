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
import os
import re
import sys

from colour import Color

import numpy as np
import pandas as pd

from generollup import __version__

_SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
_GENE_SYMBOL = "SNPEFF_TOP_EFFECT_GENE_SYMBOL"

_LOCI_COUNT = "distinct loci"
_SAMPLE_COUNT = "total impacted samples"
_MUTATION_COUNT = "total mutations"

_DBNSFP_COLUMN = "dbNSFP_rollup_damaging"
_SNPEFF_COLUMN = "SNPEFF_TOP_EFFECT_IMPACT"
_DESIRED_ANNOTATIONS = set()

_GENE_SYMBOL_OUTPUT_NAME = "gene symbol"

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

class dbNSFP(object):
    #pylint: disable=invalid-name, too-few-public-methods
    def __init__(self, format_rules, args):
        self.name = "dbNSFP"
        self.args = args
        self.column_label = "damaging votes"
        self.damaging_column = _DBNSFP_COLUMN
        self.damaging_rank_column = "overall damaging rank"
        self.damaging_total = "damaging total"
        self.format_rules = format_rules

    def summarize(self, initial_df):
        def to_rounded_string(cell_value):
            try:
                return str(int(cell_value))
            except ValueError:
                return ""

        condensed_df = self._remove_unnecessary_columns(initial_df)
        validated_df = self._remove_invalid_rows(condensed_df)

        damaging_votes_df = self._build_gene_sample_damaging_votes(validated_df)
        damaging_ranks_df = self._build_damaging_ranks(damaging_votes_df)

        frames = [damaging_ranks_df, damaging_votes_df]
        data_df = pd.concat(frames, axis=1)

        data_df = data_df.applymap(to_rounded_string)
        data_df = data_df.replace(to_replace=0,value="")

        style_dfs = []
        for format_rule in self.format_rules:
            style_dfs.append(format_rule.format(data_df))

        return data_df, style_dfs

    def _remove_unnecessary_columns(self, initial_df):
        sample_cols = initial_df.filter(regex=self.args.sample_column_regex)
        required_columns = list(sample_cols.columns)
        required_columns.extend([self.args.gene_column_name,
                                 self.damaging_column])

        condensed_df = pd.DataFrame()
        for col in initial_df.columns:
            if col in required_columns:
                condensed_df[col] = initial_df[col]

        return condensed_df

    def _remove_invalid_rows(self, initial_df):
        initial_df.fillna(".", inplace=True)
        initial_df = initial_df[initial_df[self.damaging_column] != "."]

        #temporarily reset warning thresholds to prevent spurious warning in pandas
        default = pd.options.mode.chained_assignment
        try:
            pd.options.mode.chained_assignment = None
            initial_df[self.damaging_column] = initial_df.ix[:,self.damaging_column].apply(int)
            initial_df = initial_df[initial_df[self.damaging_column] > 0]
            initial_df[self.damaging_column] = initial_df[self.damaging_column].apply(str)
        finally:
            pd.options.mode.chained_assignment = default

        return initial_df

    def _build_gene_sample_damaging_votes(self, initial_df):
        sample_cols = initial_df.filter(regex=self.args.sample_column_regex).columns.values
        initial_df = initial_df.applymap(str)

        for sample in sample_cols:
            initial_df[sample] = initial_df[sample].replace("0", ".")

            initial_df.fillna(".", inplace=True)
            initial_df[sample][initial_df[sample] != "."] = initial_df[self.damaging_column]

            split_sample = sample.split("|")
            samp_suffix = "|".join(split_sample[1:])
            sample_name = "|".join([self.name,self.column_label,samp_suffix])
            initial_df[sample] = initial_df[sample].apply(lambda x: x if x != "." else np.nan)
            initial_df[sample_name] = initial_df[sample].apply(float)
            del initial_df[sample]

        damaging_votes_df = initial_df.convert_objects(convert_numeric=True)
        damaging_votes_df = damaging_votes_df.applymap(lambda x: x if x != "." else 0)
        damaging_votes_df = damaging_votes_df.groupby(self.args.gene_column_name).sum()
        damaging_votes_df = damaging_votes_df.applymap(lambda x: x if x != 0 else ".")
        del damaging_votes_df[self.damaging_column]

        return damaging_votes_df

    def _build_damaging_ranks(self,damaging_votes_df):
        damaging_rank_name = "|".join([self.name, self.damaging_rank_column])
        total_damaging_name = "|".join([self.name, self.damaging_total])

        damaging_ranks_df = pd.DataFrame()
        damaging_ranks_df[total_damaging_name] = damaging_votes_df.sum(axis=1, skipna=True)
        damaging_ranks_df[total_damaging_name].fillna(0, inplace=True)

        damaging_ranks_df[damaging_rank_name] = damaging_ranks_df[total_damaging_name].rank(ascending=False,
                                                                                            method="min")
        damaging_ranks_df[damaging_rank_name] = damaging_ranks_df[damaging_rank_name].apply(np.floor)
        damaging_ranks_df = damaging_ranks_df.ix[:,[damaging_rank_name, total_damaging_name]]

        return damaging_ranks_df

class SnpEff(object):
    #pylint: disable=invalid-name, too-few-public-methods
    _RANK_SCORES = {"HIGH": 10**9,
                    "MODERATE": 10**6,
                    "LOW": 10**3,
                    "MODIFIER": 10**0}
    _RANK_ABBREVS = {"HIGH": "h",
                     "MODERATE": "m",
                     "LOW": "l",
                     "MODIFIER": "x",
                     ".": ""}

    def __init__(self, format_rules, args):
        self.name = "SnpEff"
        self.args = args
        self.column_label = "impact"
        self.impact_column = _SNPEFF_COLUMN
        self.impact_rank_column = "overall impact rank"
        self.impact_score_column = "overall impact score"
        self.impact_category = "impact category"
        self.format_rules = format_rules

    def summarize(self, initial_df):
        condensed_df = self._remove_unnecessary_columns(initial_df)
        validated_df = self._remove_invalid_rows(condensed_df)

        scored_df = self._calculate_score(validated_df)
        category_df = self._get_impact_category_counts(scored_df)
        ranked_df = self._calculate_rank(category_df)
        data_df = self._change_col_order(ranked_df)

        style_dfs = []
        for format_rule in self.format_rules:
            style_dfs.append(format_rule.format(data_df))

        return data_df, style_dfs

    def _remove_unnecessary_columns(self, initial_df):
        sample_cols = initial_df.filter(regex=self.args.sample_column_regex)
        required_columns = list(sample_cols.columns)
        required_columns.extend([self.args.gene_column_name,
                                 self.impact_column])

        condensed_df = pd.DataFrame()
        for col in initial_df.columns:
            if col in required_columns:
                condensed_df[col] = initial_df[col]

        return condensed_df

    def _remove_invalid_rows(self, initial_df):
        initial_df.fillna(".", inplace=True)
        initial_df = initial_df[initial_df[self.impact_column] != "."]

        return initial_df

    def _calculate_score(self, initial_df):
        sample_regex = self.args.sample_column_regex
        sample_cols = initial_df.filter(regex=sample_regex).columns.values
        initial_df = initial_df.applymap(str)

        scored_df = pd.DataFrame()
        gene_column = self.args.gene_column_name
        scored_df[gene_column] = initial_df[gene_column]

        for sample in sample_cols:
            #set sample columns equal to impact column value
            initial_df[sample][initial_df[sample] == '0'] = "."

            initial_df.fillna(".", inplace=True)
            initial_df[sample][initial_df[sample] != "."] = initial_df[self.impact_column]

            scored_df[sample + "_SCORE"] = initial_df[sample].map(SnpEff._RANK_SCORES)
            initial_df[sample] = initial_df[sample].map(SnpEff._RANK_ABBREVS)

        initial_df.fillna("", inplace=True)
        scored_df.fillna(0, inplace=True)

        score = scored_df.groupby(self.args.gene_column_name).sum().apply(sum,1)
        score = score.apply(int)
        grouped_df = initial_df.groupby(self.args.gene_column_name).sum()
        grouped_df = grouped_df[grouped_df.index != "."]
        grouped_df = grouped_df.applymap(lambda x: "".join(sorted(x)))

        if self.impact_column in grouped_df.columns.values:
            del grouped_df[self.impact_column]

        impact_score_name = "|".join([self.name, self.impact_score_column])
        grouped_df[impact_score_name] = score
        grouped_df = grouped_df.applymap(str)

        return grouped_df

    def _get_impact_category_counts(self, initial_df):
        def _count_category(abbrev):
            sample_regex = self.args.sample_column_regex
            sample_cols = initial_df.filter(regex=sample_regex).columns.values
            temp_df = pd.DataFrame()
            for sample_col in sample_cols:
                temp_df[sample_col] = initial_df[sample_col].apply(lambda x: len([i for i in x if i == abbrev]))
            return temp_df.apply(sum, 1)

        impact_category_name = "|".join([self.name, self.impact_category])
        category_df = initial_df.copy()
        category_df[impact_category_name + "|HIGH"] = _count_category("h")
        category_df[impact_category_name + "|MODERATE"] = _count_category("m")
        category_df[impact_category_name + "|LOW"] = _count_category("l")
        category_df[impact_category_name + "|MODIFIER"] = _count_category("x")

        return category_df

    def _rename_sample_columns(self, initial_df):
        sample_regex = self.args.sample_column_regex
        sample_cols = initial_df.filter(regex=sample_regex).columns.values
        for sample in sample_cols:
            split_sample = sample.split("|")
            samp_suffix = "|".join(split_sample[1:])
            sample_name = "|".join([self.name,
                                    self.column_label,
                                    samp_suffix])
            initial_df[sample_name] = initial_df[sample]
            del initial_df[sample]

        return initial_df

    def _calculate_rank(self, initial_df):
        impact_rank_name = "|".join([self.name, self.impact_rank_column])
        impact_score_name = "|".join([self.name, self.impact_score_column])

        category_df = initial_df
        category_df[impact_score_name] = category_df[impact_score_name].apply(int)
        category_df[impact_rank_name] = category_df[impact_score_name].rank(ascending=False,
                                                                            method="min")
        category_df[impact_rank_name] = category_df[impact_rank_name].apply(int)
        category_df = category_df.applymap(str)

        return self._rename_sample_columns(category_df)

    def _change_col_order(self, initial_df):
        impact_rank_name = "|".join([self.name, self.impact_rank_column])
        impact_score_name = "|".join([self.name, self.impact_score_column])

        rank = []
        score = []
        category = []
        impact = []

        for col in initial_df.columns.values:
            if col == impact_rank_name:
                rank.append(col)
            elif col == impact_score_name:
                score.append(col)
            else:
                regex = r"^{}\|{}\|.*".format(self.name, self.column_label)
                if re.match(regex, col):
                    impact.append(col)
                else:
                    category.append(col)

        ordered_names = rank + score + category + impact
        initial_df = initial_df[ordered_names]

        return initial_df


class SummaryColumns(object):
    def __init__(self, args):
        self.name = "Summary Columns"
        self.args = args

    def summarize(self, initial_df):
        sample_df = initial_df.filter(regex=self.args.sample_column_regex)
        gene_column = self.args.gene_column_name
        sample_df[gene_column] = initial_df[gene_column]

        data_df = pd.DataFrame()
        data_df[_SAMPLE_COUNT] = self.calculate_total_samples(sample_df)
        data_df[_LOCI_COUNT] = self.calculate_total_loci(sample_df)
        data_df[_MUTATION_COUNT] = self.calculate_total_mutations(sample_df)
        data_df = data_df.applymap(str)

        style_dfs = []
        style_df = data_df.copy()
        style_dfs.append(style_df.applymap(lambda x: ""))

        return data_df, style_dfs

    def calculate_total_samples(self, sample_df):
        sample_df[sample_df == '.'] = np.nan
        sample_df[sample_df == '0'] = np.nan
        sample_df = sample_df.groupby(self.args.gene_column_name).count()

        sample_df[sample_df > 0] = 1

        return sample_df.apply(sum, 1)

    def calculate_total_loci(self, sample_df):
        return sample_df.groupby(self.args.gene_column_name).size()

    def calculate_total_mutations(self, sample_df):
        sample_df[sample_df == '.'] = np.nan
        sample_df[sample_df == '0'] = np.nan
        grouped = sample_df.groupby(self.args.gene_column_name)
        return grouped.count().apply(sum, 1)


class SnpEffFormatRule(object):
    #pylint: disable=too-few-public-methods
    _COLOR_MAP = {"h": "#003366",
                  "m": "#3377B9",
                  "l": "#91C4E8",
                  "x": "#CCCCCC"}
    _FONT_COLOR = "#000000"

    def __init__(self):
        self.rank_column = "SnpEff|overall impact rank"
        self.range_rule = list(Color("white").range_to(Color("blue"), 101))

    def format(self, data_df):
        def _determine_format(cell_value):
            for letter in SnpEffFormatRule._COLOR_MAP.keys():
                if letter in cell_value:
                    return letter
            return ""

        format_df = data_df.applymap(str)
        format_df = format_df.applymap(_determine_format)
        format_df = format_df.applymap(self._style)

        return format_df

    @staticmethod
    def _style(cell_value):
        if len(cell_value) < 1:
            return np.nan

        else:
            styled_cell = {"font_size": "4",
                           "bg_color": SnpEffFormatRule\
                                               ._COLOR_MAP[cell_value],
                           "font_color": SnpEffFormatRule\
                                         ._COLOR_MAP[cell_value]}
            return styled_cell

#     def _style_score(self, cell_value):
#         if np.isnan(cell_value):
#             return np.nan
#         else:
#             color = self.range_rule[cell_value]
#             styled_cell = {"font_size": "12",
#                            "bg_color": str(color),
#                            "font_color": SnpEffFormatRule._FONT_COLOR}
#             return styled_cell


class dbNSFPFormatRule(object):
    #pylint: disable=invalid-name,too-few-public-methods

    _FONT_COLOR = "#000000"

    def __init__(self):
        self.rank_column = "dbNSFP|overall damaging rank"
        self.total_column = "dbNSFP|damaging total"
        self.range_rule = list(Color("white").range_to(Color("orange"), 101))

    def format(self, data_df):
        def _determine_format(cell_value):
            normalized = 100*(cell_value-min_value)/(max_value-min_value)
            if np.isnan(normalized):
                return cell_value
            else:
                return int(normalized)

        format_df = data_df.drop(self.rank_column, 1)
        format_df = format_df.drop(self.total_column, 1)

        format_df = format_df.replace("", np.nan)
        format_df = format_df.convert_objects(convert_numeric=True)

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
            styled_cell = {"font_size": "12",
                           "bg_color": str(color),
                           "font_color": dbNSFPFormatRule._FONT_COLOR}
            return styled_cell


class RankFormatRule(object):
    #pylint: disable=invalid-name, too-few-public-methods

    _FONT_COLOR = "#000000"

    def __init__(self):
        self.range_rule = list(Color("red").range_to(Color("white"), 101))

    def format(self, data_df):
        def _determine_format(cell_value):
            if len(str(cell_value)) > 0:
                cell_value = int(cell_value)
                normalized = 100*(cell_value-min_value)/(max_value-min_value)
                return str(int(normalized))
            else:
                return ""

        format_df = data_df.copy()

        rank_columns = {}

        for col in format_df.columns.values:
            if re.search(dbNSFP("","").damaging_rank_column, col) or \
                re.search(SnpEff("","").impact_rank_column, col):
                min_value = int(format_df[col].min(skipna=True))
                max_value = int(format_df[col].max(skipna=True))

                if max_value - min_value > 0:
                    format_df.fillna("", inplace=True)
                    format_df[col] = format_df[col].apply(_determine_format)

                    rank_columns[col] = format_df[col]

                else:
                    format_df[col] = format_df[col].apply(lambda x: "")

            else:
                format_df[col] = format_df[col].apply(lambda x: "")

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
        self.range_rule = list(Color("#E27171").range_to(Color("white"), value))

    def _style(self, cell_value):
        if len(str(cell_value)) > 0:
            cell_value = int(cell_value)
            color = self.range_rule[cell_value]
            styled_cell = {"font_size": "12",
                           "bg_color": str(color),
                           "font_color": RankFormatRule._FONT_COLOR}
            return styled_cell
        return np.nan


def _create_df(input_file, args):
    initial_df = pd.read_csv(input_file, sep='\t', header=0, dtype='str')
    _validate_df(initial_df, args)
    initial_df = initial_df[pd.notnull(initial_df[_GENE_SYMBOL])]

    return initial_df

def _validate_df(initial_df, args):
    header = set(initial_df.columns.values)
    missing_columns = []

    if _DBNSFP_COLUMN in header:
        _DESIRED_ANNOTATIONS.add(_DBNSFP_COLUMN)
    if _SNPEFF_COLUMN in header:
        _DESIRED_ANNOTATIONS.add(_SNPEFF_COLUMN)
    elif _DBNSFP_COLUMN not in header and _SNPEFF_COLUMN not in header:
        missing_columns.extend([_DBNSFP_COLUMN, _SNPEFF_COLUMN])

    if args.gene_column_name not in header:
        missing_columns.append(args.gene_column_name)

    msg = ("Input file is missing required headers ({}). "
           "Review input and try again.").format(missing_columns)
    if missing_columns:
        raise UsageError(msg)

    sample_column = 0
    for column in header:
        if re.search(args.sample_column_regex, column):
            sample_column = 1
            break

    msg = ("Cannot determine samples from input file with supplied regex. "
           "Review input and try again.")
    if not sample_column:
        raise UsageError(msg)

def _combine_style_dfs(dfs):
    df1, df2 = dfs
    original_columns = set(df1.columns.values).union(set(df2.columns.values))

    cleaned_df1 = df1.dropna(1, how='all')
    cleaned_df2 = df2.dropna(1, how='all')

    combined_df = cleaned_df1.join(cleaned_df2, how='outer')
    for column in original_columns:
        if column not in list(combined_df.columns.values):
            combined_df[column] = pd.Series()

    combined_df = combined_df.fillna("")

    return combined_df

def _combine_dfs(dfs):
    summary_df = dfs.values()[0]
    annotation_dfs = dfs.values()[1:]

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

    combined_df = combined_df.fillna("")
    combined_df.index.names = [_GENE_SYMBOL_OUTPUT_NAME]

    combined_df = combined_df[combined_df.index != "."]

    return combined_df

def _header_formats():
    #pylint: disable=duplicate-key
    dbnsfp = dbNSFP("","")
    snpeff = SnpEff("","")
    header_formats = {_GENE_SYMBOL_OUTPUT_NAME: {"align": "center",
                                                 "bold": True,
                                                 "border": True,
                                                 "align": "center"},
                      dbnsfp.damaging_rank_column: {"bg_color": "#F8A049",
                                                    "font_color": "white",
                                                    "border": True,
                                                    "rotation": 90,
                                                    "align": "center"},
                      dbnsfp.damaging_total: {"bg_color": "#FBC692",
                                              "border": True,
                                              "rotation": 90,
                                              "align": "center"},
                      dbnsfp.column_label: {"bg_color": "#FBC692",
                                            "border": True,
                                            "rotation": 90,
                                            "align": "center"},
                      snpeff.impact_rank_column: {"bg_color": "#4775A3",
                                                  "font_color": "white",
                                                  "border": True,
                                                  "rotation": 90,
                                                  "align": "center"},
                      snpeff.impact_score_column: {"bg_color": "#6C91B5",
                                                   "font_color": "white",
                                                   "border": True,
                                                   "rotation": 90,
                                                   "align": "center"},
                      snpeff.impact_category: {"bg_color": "#6C91B5",
                                               "font_color": "white",
                                               "border": True,
                                               "rotation": 90,
                                               "align": "center"},
                      snpeff.column_label + r"\|": {"bg_color": "#99C1D7",
                                                    "border": True,
                                                    "rotation": 90,
                                                    "align": "center"},
                      _LOCI_COUNT: {"bg_color": "#F7F4FF",
                                    "border": True,
                                    "rotation": 90,
                                    "align": "center"},
                      _SAMPLE_COUNT: {"bg_color": "#F7F4FF",
                                      "border": True,
                                      "rotation": 90,
                                      "align": "center"},
                      _MUTATION_COUNT: {"bg_color": "#F7F4FF",
                                        "border": True,
                                        "rotation": 90,
                                        "align": "center"}}
    return header_formats

def _translate_to_excel(data_df, style_df, writer):
    gene_symbols = data_df[_GENE_SYMBOL_OUTPUT_NAME].values
    data_df = data_df.convert_objects(convert_numeric=True)
    data_df[_GENE_SYMBOL_OUTPUT_NAME] = gene_symbols

    worksheet_name = "gene_rollup"
    data_df.to_excel(writer, sheet_name=worksheet_name, index=False)

    workbook = writer.book
    workbook.strings_to_numbers = True
    worksheet = writer.sheets[worksheet_name]
    worksheet.strings_to_numbers = True

    worksheet.set_column(0, 0, 20)
    worksheet.set_column(1, len(list(data_df.columns.values)), 5)

    for i, column in enumerate(data_df.columns.values):
        for key in _header_formats():
            if re.search(dbNSFP("","").damaging_total, column) or \
                re.search(SnpEff("","").impact_score_column, column):
                worksheet.set_column(i, i, None, None, {"hidden": 1})

            if re.search(key, column):
                cell_format = workbook.add_format(_header_formats()[key])
                worksheet.write(0, i, column, cell_format)
                break

    print "{} | Writing data to an Excel file".format(_now())
    for i, (row, dummy) in enumerate(data_df.iterrows()):
        for j, (column, dummy) in enumerate(data_df.iteritems()):
            style = style_df.ix[row, column]
            if style:
                cell_format = workbook.add_format(style)
                worksheet.write(i + 1, j, data_df.ix[row, column], cell_format)


    writer.save()

def _sort_by_dbnsfp_rank(initial_df):
    dbnsfp = dbNSFP("","")
    snpeff = SnpEff("","")
    dbnsfp_rank_col = "|".join([dbnsfp.name, dbnsfp.damaging_rank_column])
    snpeff_rank_col = "|".join([snpeff.name, snpeff.impact_rank_column])

    sort_order = []
    try:
        initial_df[dbnsfp_rank_col] = initial_df[dbnsfp_rank_col].replace("",
                                                                          np.nan)
        initial_df[dbnsfp_rank_col] = initial_df[dbnsfp_rank_col].apply(float)
        sort_order.append(dbnsfp_rank_col)
    except KeyError:
        pass

    try:
        initial_df[snpeff_rank_col] = initial_df[snpeff_rank_col].replace("",
                                                                          np.nan)
        initial_df[snpeff_rank_col] = initial_df[snpeff_rank_col].apply(float)
        sort_order.append(snpeff_rank_col)
    except KeyError:
        pass

    sort_order.append(_GENE_SYMBOL_OUTPUT_NAME)
    sorted_df = initial_df.sort(columns = sort_order)

    sorted_df.fillna("", inplace=True)

    return sorted_df

def _reset_style_gene_values(combined_style_df):
    original_values = combined_style_df[_GENE_SYMBOL_OUTPUT_NAME]
    new_values = [{} for dummy in original_values]
    combined_style_df[_GENE_SYMBOL_OUTPUT_NAME] = new_values

    return combined_style_df

def _rollup(args):
    initial_df = _create_df(args.input_file, args)
    print "{} | Starting Gene Rollup".format(_now())

    annotations = [SummaryColumns(args)]
    if _DBNSFP_COLUMN in _DESIRED_ANNOTATIONS:
        annotations.append(dbNSFP([dbNSFPFormatRule(), RankFormatRule()], args))
    if _SNPEFF_COLUMN in _DESIRED_ANNOTATIONS:
        annotations.append(SnpEff([SnpEffFormatRule(), RankFormatRule()], args))

    annotation_dfs = OrderedDict()
    all_style_dfs = OrderedDict()

    for annotation in annotations:
        print "{} | Generating {} rollup information".format(_now(),
                                                             annotation.name)
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

    print "{} | Writing to output file [{}]".format(_now(), args.output_file)

    if args.tsv:
        sorted_df.to_csv(args.output_file, sep="\t", index=False)
    else:
        try:
            writer = pd.ExcelWriter(args.output_file, engine="xlsxwriter")
            _translate_to_excel(sorted_df, altered_style_df, writer)

        except ValueError:
            msg = ("Unable to write [{}] to an Excel file. Review inputs and "
                   "try again.").format(args.output_file)
            raise RollupException(msg)

    print "{} | Wrote to [{}]".format(_now(),
                                      args.output_file)
    print "{} | Done.".format(_now())

def _create_output_dir(output_file):
    try:
        _makepath(output_file)
    except OSError:
        parent_dir = os.path.dirname(output_file)
        raise UsageError(("GeneRollup cannot write to output directory "
                          "[{}]. Review inputs and try again.")\
                          .format(parent_dir))

def _makepath(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def _now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def _add_arg_parse(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_file", help=("A TSV file of variants x samples"))
    parser.add_argument("output_file",
                        help=("By default, an XLSX file of genes x samples. "
                              "Can be a TSV file with --tsv"))
    parser.add_argument("--sample_column_regex",
                        default=_SAMPLENAME_REGEX,
                        help=("Regex used to define the sample columns in the "
                              "input file"))
    parser.add_argument("--gene_column_name",
                        default=_GENE_SYMBOL,
                        help="Name of gene symbol column in the output file")
    parser.add_argument("--tsv",
                        action="store_true",
                        default=False,
                        help="Write to a tsv file rather than an xlsx file")
    parser.add_argument("--version",
                        "-V",
                        action="version",
                        version=__version__)
    return parser.parse_args(args)

def main():
    args = _add_arg_parse(sys.argv[1:])
#     TODO: (jebene) make this work
#     _create_output_dir(args.output_file)

    try:
        _rollup(args)
    except UsageError as error:
        print "{} | Usage Error: {}".format(_now(), error.message)

if __name__ == "__main__":
    main()
