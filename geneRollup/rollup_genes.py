#!/usr/bin/env python
#pylint:disable=line-too-long, no-member
import argparse
from collections import OrderedDict
import re
import sys

from colour import Color

import pandas as pd
import numpy as np


_SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
_GENE_SYMBOL = "SNPEFF_TOP_EFFECT_GENE_SYMBOL"
_XLSX = True

_LOCI_COUNT = "distinct loci"
_SAMPLE_COUNT = "total impacted samples"
_MUTATION_COUNT = "total mutations"

_REQUIRED_COLUMNS = set(["dbNSFP_rollup_damaging",
                         "SNPEFF_TOP_EFFECT_IMPACT"])

_GENE_SYMBOL_OUTPUT_NAME = "gene symbol"

class dbNSFP(object):
    #pylint: disable=invalid-name, too-few-public-methods
    def __init__(self, format_rules):
        self.name = "dbNSFP"
        self.column_label = "damaging votes"
        self.damaging_column = "dbNSFP_rollup_damaging"
        self.damaging_rank_column = "dbNSFP|overall damaging rank"
        self.damaging_total = "dbNSFP|damaging total"
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
        sample_cols = initial_df.filter(regex=_SAMPLENAME_REGEX)
        required_columns = list(sample_cols.columns)
        required_columns.extend([_GENE_SYMBOL, self.damaging_column])

        condensed_df = pd.DataFrame()
        for col in initial_df.columns:
            if col in required_columns:
                condensed_df[col] = initial_df[col]

        return condensed_df

    def _remove_invalid_rows(self, initial_df):
        initial_df.fillna(".", inplace=True)
        initial_df = initial_df[initial_df[self.damaging_column] != "."]

        initial_df[self.damaging_column] = initial_df[self.damaging_column].apply(int)
        initial_df = initial_df[initial_df[self.damaging_column] > 0]

        initial_df[self.damaging_column] = initial_df[self.damaging_column].apply(str)

        return initial_df

    def _build_gene_sample_damaging_votes(self, initial_df):
        sample_cols = initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
        initial_df = initial_df.applymap(str)

        for sample in sample_cols:
            initial_df[sample][initial_df[sample] == '0'] = "."

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
        damaging_votes_df = damaging_votes_df.groupby(_GENE_SYMBOL).sum()
        damaging_votes_df = damaging_votes_df.applymap(lambda x: x if x != 0 else ".")
        del damaging_votes_df[self.damaging_column]

        return damaging_votes_df

    def _build_damaging_ranks(self,damaging_votes_df):
        damaging_votes_df = damaging_votes_df.filter(regex="dbNSFP.*TUMOR")
        damaging_ranks_df = pd.DataFrame()
        damaging_ranks_df[self.damaging_total] = damaging_votes_df.sum(axis=1, skipna=True)
        damaging_ranks_df[self.damaging_total].fillna(0, inplace=True)

        damaging_ranks_df[self.damaging_rank_column] = damaging_ranks_df[self.damaging_total].rank(ascending=False)
        damaging_ranks_df[self.damaging_rank_column] = damaging_ranks_df[self.damaging_rank_column].apply(np.floor)
        damaging_ranks_df = damaging_ranks_df.ix[:,[self.damaging_rank_column,self.damaging_total]]

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

    def __init__(self, format_rules):
        self.name = "SnpEff"
        self.column_label = "impact"
        self.impact_column = "SNPEFF_TOP_EFFECT_IMPACT"
        self.impact_rank_column = "SnpEff|overall impact rank"
        self.impact_score_column = "SnpEff|overall impact score"
        self.impact_category = "SnpEff|impact category|{}"
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
        sample_cols = initial_df.filter(regex=_SAMPLENAME_REGEX)
        required_columns = list(sample_cols.columns)
        required_columns.extend([_GENE_SYMBOL, self.impact_column])

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
    #pylint: disable=line-too-long
        sample_cols = initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
        initial_df = initial_df.applymap(str)

        scored_df = pd.DataFrame()
        scored_df[_GENE_SYMBOL] = initial_df[_GENE_SYMBOL]

        for sample in sample_cols:
            #set sample columns equal to impact column value
            initial_df[sample][initial_df[sample] == '0'] = "."

            initial_df.fillna(".", inplace=True)
            initial_df[sample][initial_df[sample] != "."] = initial_df[self.impact_column]

            scored_df[sample + "_SCORE"] = initial_df[sample].map(SnpEff._RANK_SCORES)
            initial_df[sample] = initial_df[sample].map(SnpEff._RANK_ABBREVS)

        initial_df.fillna("", inplace=True)
        scored_df.fillna(0, inplace=True)

        score = scored_df.groupby(_GENE_SYMBOL).sum().apply(sum, 1)
        score = score.apply(int)
        grouped_df = initial_df.groupby(_GENE_SYMBOL).sum()
        grouped_df = grouped_df[grouped_df.index != "."]
        grouped_df = grouped_df.applymap(lambda x: "".join(sorted(x)))

        if self.impact_column in grouped_df.columns.values:
            del grouped_df[self.impact_column]

        grouped_df[self.impact_score_column] = score
        grouped_df = grouped_df.applymap(str)

        return grouped_df

    def _get_impact_category_counts(self, initial_df):
        #pylint: disable=line-too-long
        def _count_category(abbrev):
            sample_cols = initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
            temp_df = pd.DataFrame()
            for sample_col in sample_cols:
                temp_df[sample_col] = initial_df[sample_col].apply(lambda x: len([i for i in x if i == abbrev]))
            return temp_df.apply(sum, 1)

        category_df = initial_df.copy()
        category_df[self.impact_category.format("HIGH")] = _count_category("h")
        category_df[self.impact_category.format("MODERATE")] = _count_category("m")
        category_df[self.impact_category.format("LOW")] = _count_category("l")
        category_df[self.impact_category.format("MODIFIER")] = _count_category("x")

        return category_df

    def _rename_sample_columns(self, initial_df):
        sample_cols = initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
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
        #pylint: disable=line-too-long
        category_df = initial_df
        category_df[self.impact_score_column] = category_df[self.impact_score_column].apply(int)
        category_df[self.impact_rank_column] = category_df[self.impact_score_column].rank(ascending=0, method="min")
        category_df[self.impact_rank_column] = category_df[self.impact_rank_column].apply(int)
        category_df = category_df.applymap(str)
        return self._rename_sample_columns(category_df)

    def _change_col_order(self, initial_df):
        rank = []
        score = []
        category = []
        impact = []
        for col in initial_df.columns.values:
            if col == self.impact_rank_column:
                rank.append(col)
            elif col == self.impact_score_column:
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
    def __init__(self):
        self.name = "Summary Columns"

    def summarize(self, initial_df):
        sample_df = initial_df.filter(regex=_SAMPLENAME_REGEX)
        sample_df[_GENE_SYMBOL] = initial_df[_GENE_SYMBOL]

        data_df = pd.DataFrame()
        data_df[_SAMPLE_COUNT] = self.calculate_total_samples(sample_df)
        data_df[_LOCI_COUNT] = self.calculate_total_loci(sample_df)
        data_df[_MUTATION_COUNT] = self.calculate_total_mutations(sample_df)
        data_df = data_df.applymap(str)

        style_dfs = []
        style_df = data_df.copy()
        style_dfs.append(style_df.applymap(lambda x: ""))

        return data_df, style_dfs

    @staticmethod
    def calculate_total_samples(sample_df):
        sample_df[sample_df == '.'] = np.nan
        sample_df[sample_df == '0'] = np.nan
        sample_df = sample_df.groupby(_GENE_SYMBOL).count()
        sample_df[sample_df > 0] = 1

        return sample_df.apply(sum, 1)

    @staticmethod
    def calculate_total_loci(sample_df):
        return sample_df.groupby(_GENE_SYMBOL).size()

    @staticmethod
    def calculate_total_mutations(sample_df):
        sample_df[sample_df == '.'] = np.nan
        sample_df[sample_df == '0'] = np.nan
        grouped = sample_df.groupby(_GENE_SYMBOL)
        return grouped.count().apply(sum, 1)


class SnpEffFormatRule(object):
    #pylint: disable=too-few-public-methods
    _COLOR_MAP = {"h": "#003366",
                  "m": "#6699FF",
                  "l": "#99CCFF",
                  "x": "#CCCCCC"}

    def format(self, data_df):
        def _determine_format(cell_value):
            for letter in SnpEffFormatRule._COLOR_MAP.keys():
                if letter in cell_value:
                    return letter
            return ""

        format_df = data_df.applymap(str)
        format_df = format_df.applymap(_determine_format)

        for column, dummy in format_df.iteritems():
            for row, dummy in format_df.iterrows():
                styled_cell = self._style(format_df.ix[row, column])
                format_df.ix[row, column] = styled_cell

        return format_df

    @staticmethod
    def _style(cell_value):
        if len(cell_value) > 0:
            styled_cell = {"font_size": "4",
                           "bg_color": SnpEffFormatRule\
                                               ._COLOR_MAP[cell_value],
                           "font_color": SnpEffFormatRule\
                                         ._COLOR_MAP[cell_value]}
            return styled_cell
        return ""


class dbNSFPFormatRule(object):
    #pylint: disable=invalid-name,too-few-public-methods

    _RANGE_RULE = list(Color("white").range_to(Color("orange"), 101))
    _FONT_COLOR = "#000000"

    def __init__(self):
        self.rank_column = "dbNSFP|overall damaging rank"

    def format(self, data_df):
        def _determine_format(cell_value):
            normalized = 100*(cell_value-min_value)/(max_value-min_value)
            if np.isnan(normalized):
                return cell_value
            else:
                return int(normalized)

        #TODO: (jebene) make "dbNSFP|overall damaging rank" a constant -pass in?
        format_df = data_df.drop(self.rank_column, 1)
        format_df = format_df.replace("", np.nan)
        format_df = format_df.convert_objects(convert_numeric=True)

        min_value = min(format_df.min(skipna=True))
        max_value = max(format_df.max(skipna=True))

#         format_df.fillna("", inplace=True)
        format_df = format_df.applymap(_determine_format)

        format_df[self.rank_column] = pd.Series()
#         format_df.fillna("", inplace=True)

        format_df = format_df.applymap(self._style)
#         for column, dummy in format_df.iteritems():
#             for row, dummy in format_df.iterrows():
#                 styled_cell = self._style(format_df.ix[row, column])
#                 format_df.ix[row, column] = styled_cell
        return format_df

    @staticmethod
    def _style(cell_value):
        if np.isnan(cell_value):
            return {}
        else:
            cell_value = int(cell_value)
            color = dbNSFPFormatRule._RANGE_RULE[cell_value]
            styled_cell = {"font_size": "12",
                           "bg_color": str(color),
                           "font_color": dbNSFPFormatRule._FONT_COLOR}
            return styled_cell


#TODO: (jebene) hookup rankformatrule
class RankFormatRule(object):
    #pylint: disable=invalid-name, too-few-public-methods
    _RANGE_RULE = list(Color("red").range_to(Color("white"), 101))
    _FONT_COLOR = "#000000"

    def __init__(self):
        #TODO: (jebene) modify regex storage so that annotations and format rule refer to same format
        self.rank_regex = r".*\|overall .* rank"

    def format(self, data_df):
        def _determine_format(cell_value):
            if len(str(cell_value)) > 0:
                normalized = 100*(cell_value-min_value)/(max_value-min_value)
                return str(int(normalized))
            else:
                return ""

        format_df = data_df.copy()
        for col in format_df.columns.values:
            if re.search(self.rank_regex, col):
                min_value = int(format_df[col].min(skipna=True))
                max_value = int(format_df[col].max(skipna=True))

                format_df.fillna("", inplace=True)
                format_df[col] = format_df[col].apply(_determine_format)
            else:
                format_df[col] = format_df[col].apply(lambda x: "")

        for column, dummy in format_df.iteritems():
            for row, dummy in format_df.iterrows():
                styled_cell = self._style(format_df.ix[row, column])
                format_df.ix[row, column] = styled_cell

        return format_df

    @staticmethod
    def _style(cell_value):
        if len(str(cell_value)) > 0:
            cell_value = int(cell_value)
            color = RankFormatRule._RANGE_RULE[cell_value]
            styled_cell = {"font_size": "12",
                           "bg_color": str(color),
                           "font_color": RankFormatRule._FONT_COLOR}
            return styled_cell
        return ""

def _create_df(input_file):
    initial_df = pd.read_csv(input_file, sep='\t', header=False, dtype='str')
    _validate_df(initial_df)

    return initial_df

def _validate_df(initial_df):
    header = set(initial_df.columns.values)
    _REQUIRED_COLUMNS.add(_GENE_SYMBOL)
    missing_columns = _REQUIRED_COLUMNS.difference(header)
    msg = ("Input file is missing required headers ({}). "
           "Review input and try again.").format(missing_columns)
    if missing_columns:
        raise BaseException(msg)

    msg = ("Cannot determine samples from input file with supplied regex. "
           "Review input and try again.")
    sample_column = 0
    for column in header:
        if re.search(_SAMPLENAME_REGEX, column):
            sample_column = 1
            break
    if not sample_column:
        raise BaseException(msg)

#TODO: (jebene) - edit this so that only columns with data are joined
def _combine_style_dfs(dfs):
    df1, df2 = dfs
    combined_df = df1.join(df2, how='outer')
    combined_df = combined_df.fillna("")
    return combined_df

def _combine_dfs(dfs):
    summary_df, dbnsfp_df, snpeff_df = dfs.values()

    dbnsfp_genes = set(dbnsfp_df.index.values)
    snpeff_genes = set(snpeff_df.index.values)
    summary_genes = set(summary_df.index.values)
    valid_genes = dbnsfp_genes.intersection(snpeff_genes)
    invalid_genes = summary_genes.difference(valid_genes)
 
    summary_df = summary_df.drop(list(invalid_genes))

    combined_df = summary_df.join(dbnsfp_df, how='outer')
    combined_df = combined_df.join(snpeff_df, how='outer')
    combined_df = combined_df.fillna("")
    combined_df.index.names = [_GENE_SYMBOL_OUTPUT_NAME]

    combined_df = combined_df[combined_df.index != "."]

    return combined_df

def _translate_to_excel(data_df, style_df, writer):
    worksheet_name = "gene_rollup"
    data_df.to_excel(writer, sheet_name=worksheet_name, index=False)

    workbook = writer.book
    worksheet = writer.sheets[worksheet_name]

    for i, (row, dummy) in enumerate(data_df.iterrows()):
        for j, (column, dummy) in enumerate(data_df.iteritems()):
            style = style_df.ix[row, column]

            if len(style) > 0 and column != _GENE_SYMBOL_OUTPUT_NAME:
                cell_format = workbook.add_format(style)
                worksheet.write(i + 1, j, data_df.ix[row, column], cell_format)

    writer.save()

def _sort_by_dbnsfp_rank(initial_df):
    sorted_df = initial_df.sort(columns = [dbNSFP("").damaging_rank_column,
                                           SnpEff("").impact_rank_column,
                                           _GENE_SYMBOL_OUTPUT_NAME])
    return sorted_df

def _rollup(input_file, output_file):
    print "Starting Gene Rollup"

    initial_df = _create_df(input_file)
#TODO: (jebene) add RankFormatRule() to dbNSFP and SnpEff -- it's reusing data_df in summarize() I think, so look into that
    annotations = [SummaryColumns(),
                   dbNSFP([dbNSFPFormatRule()]),
                   SnpEff([SnpEffFormatRule()])]

    annotation_dfs = OrderedDict()
    all_style_dfs = OrderedDict()

    for annotation in annotations:
        print "Generating {} rollup information".format(annotation.name)
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

    if _XLSX:
        try:
            writer = pd.ExcelWriter(output_file, engine="xlsxwriter")
            _translate_to_excel(sorted_df, combined_style_df, writer)

        except ValueError:
            msg = ("Unable to write [{}] to an Excel file. Review inputs and "
                   "try again.").format(output_file)
            raise BaseException(msg)
    else:
        sorted_df.to_csv(output_file, sep="\t", index=False)

    print "Wrote to [{}]".format(output_file)
    print "Done."

def _change_global_variables(args):
    #pylint: disable=global-statement
    global _SAMPLENAME_REGEX
    if args.sample_column_regex:
        _SAMPLENAME_REGEX = args.sample_column_regex

    global _GENE_SYMBOL
    if args.gene_column_name:
        _GENE_SYMBOL = args.gene_column_name

    global _XLSX
    if args.tsv:
        _XLSX = False

def _add_arg_parse(args):
    parser = argparse.ArgumentParser()
    #pylint: disable=line-too-long
    parser.add_argument("input_file", help=("A tab-delimited file of variants x samples"))
    parser.add_argument("output_file", help=("A tab-delimited file of genes x samples"))
    parser.add_argument("--sample_column_regex", help="Regex used to define the sample columns in the input file. Default is 'JQ_SUMMARY_SOM_COUNT.*'")
    parser.add_argument("--gene_column_name", help="Name of gene symbol column in the output file. Default is 'SNPEFF_TOP_EFFECT_GENE_SYMBOL'")
    parser.add_argument("--tsv", action="store_true", help="Write to a tsv file rather than an xlsx file")

    return parser.parse_args(args)

def main():
    args = _add_arg_parse(sys.argv[1:])
    _change_global_variables(args)
    _rollup(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
