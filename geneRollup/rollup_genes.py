#!/usr/bin/env python
import argparse
from collections import OrderedDict
import re
import sys

from colour import Color

import pandas as pd


_SAMPLENAME_REGEX = ""
_GENE_SYMBOL = ""
_XLSX = False

_LOCI_COUNT = "distinct loci"
_SAMPLE_COUNT = "total impacted samples"
_MUTATION_COUNT = "total mutations"

_REQUIRED_COLUMNS = set(["dbNSFP_rollup_damaging",
                         "SNPEFF_TOP_EFFECT_IMPACT"])

_GENE_SYMBOL_OUTPUT_NAME = "gene symbol"

class dbNSFP(object):
    #pylint: disable=invalid-name, too-few-public-methods
    def __init__(self, format_rule):
        self.name = "dbNSFP"
        self.column_label = "damaging votes"
        self.damaging_column = "dbNSFP_rollup_damaging"
        self.damaging_rank_column = "dbNSFP|overall damaging rank"
        self.format_rule = format_rule

    def summarize(self, initial_df):
        condensed_df = self._remove_unnecessary_columns(initial_df)
        ranked_df = self._calculate_rank(condensed_df)
        data_df = self._change_col_order(ranked_df)
        style_df = self.format_rule.format(data_df)

        return data_df, style_df

    def _remove_unnecessary_columns(self, initial_df):
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX)
        required_columns = list(sample_cols.columns)
        required_columns.extend([_GENE_SYMBOL, self.damaging_column])

        condensed_df = pd.DataFrame()
        for col in initial_df.columns:
            if col in required_columns:
                condensed_df[col] = initial_df[col]

        return condensed_df

    def _calculate_rank(self, initial_df):
        #pylint: disable=line-too-long
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
        initial_df = initial_df.applymap(str)

        #set sample columns equal to damaging column value
        for sample in sample_cols:
            initial_df[sample][initial_df[sample] != "."] = initial_df[self.damaging_column]

        for sample in sample_cols:
            split_sample = sample.split("|")
            samp_suffix = "|".join(split_sample[1:])
            sample_name = "|".join([self.name,
                                    self.column_label,
                                    samp_suffix])
            initial_df[sample] = initial_df[sample].apply(lambda x: x if x != "." else 0)
            initial_df[sample_name] = initial_df[sample].apply(float)

            del initial_df[sample]

        ranked_df = initial_df.groupby(_GENE_SYMBOL).sum()
        ranked_df[self.damaging_rank_column] = ranked_df.apply(sum, 1)

        ranked_df.fillna("", inplace=True)

        ranked_df = ranked_df.sort(self.damaging_rank_column, ascending=0)
        ranked_df[self.damaging_rank_column] = ranked_df[self.damaging_rank_column].rank(ascending=0, method="min")

        return ranked_df

    def _change_col_order(self, initial_df):
        ordered_names = []
        rank = []
        for col in initial_df.columns.values:
            if col == self.damaging_rank_column:
                rank.append(col)
            else:
                ordered_names.append(col)

        all_names = rank + ordered_names

        initial_df = initial_df[all_names]
        return initial_df


class SnpEff(object):
    #pylint: disable=invalid-name, too-few-public-methods
    _RANK_SCORES =  {"HIGH": 100000.0,
                     "MODERATE": 1,
                     "LOW": 1/100000.0,
                     "MODIFIER": 1/10**12}
    _RANK_ABBREVS =  {"HIGH": "h", "MODERATE": "m", "LOW": "l", "MODIFIER": "x"}

    def __init__(self, format_rule):
        self.name = "SnpEff"
        self.column_label = "impact"
        self.impact_column = "SNPEFF_TOP_EFFECT_IMPACT"
        self.impact_rank_column = "SnpEff|overall impact rank"
        self.impact_score_column = "SnpEff|overall impact score"
        self.impact_category = "SnpEff|impact category|{}"
        self.format_rule = format_rule

    def summarize(self, initial_df):
        condensed_df = self._remove_unnecessary_columns(initial_df)
        scored_df = self._calculate_score(condensed_df)
        category_df = self._get_impact_category_counts(scored_df)
        ranked_df = self._calculate_rank(category_df)
        data_df = self._change_col_order(ranked_df)
        style_df = self.format_rule.format(data_df)

        return data_df, style_df

    def _remove_unnecessary_columns(self, initial_df):
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX)
        required_columns = list(sample_cols.columns)
        required_columns.extend([_GENE_SYMBOL, self.impact_column])

        condensed_df = pd.DataFrame()
        for col in initial_df.columns:
            if col in required_columns:
                condensed_df[col] = initial_df[col]

        return condensed_df

    def _calculate_score(self, initial_df):
    #pylint: disable=line-too-long
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
        initial_df = initial_df.applymap(str)

        scored_df = pd.DataFrame()
        scored_df[_GENE_SYMBOL] = initial_df[_GENE_SYMBOL]

        for sample in sample_cols:
            #set sample columns equal to impact column value
            initial_df[sample][initial_df[sample] != "."] = initial_df[self.impact_column]

            scored_df[sample + "_SCORE"] = initial_df[sample].map(SnpEff._RANK_SCORES)
            initial_df[sample] = initial_df[sample].map(SnpEff._RANK_ABBREVS)

        initial_df.fillna("", inplace=True)
        scored_df.fillna(0, inplace=True)

        score = scored_df.groupby(_GENE_SYMBOL).sum().apply(sum, 1)
        grouped_df = initial_df.groupby(_GENE_SYMBOL).sum()
        grouped_df = grouped_df[grouped_df.index != "."]
        grouped_df = grouped_df.applymap(lambda x: "".join(sorted(x)))

        if self.impact_column in grouped_df.columns.values:
            del grouped_df[self.impact_column]

        grouped_df[self.impact_score_column] = score

        return grouped_df

    def _get_impact_category_counts(self, initial_df):
        #pylint: disable=line-too-long
        def _count_category(abbrev):
            sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
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
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
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
        category_df = initial_df.sort(self.impact_score_column, ascending=0)
        category_df[self.impact_rank_column] = category_df[self.impact_score_column].rank(ascending=0, method="min")

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
        sample_df =  initial_df.filter(regex=_SAMPLENAME_REGEX)
        sample_df[_GENE_SYMBOL] = initial_df[_GENE_SYMBOL]

        data_df = pd.DataFrame()
        data_df[_SAMPLE_COUNT] = self.calculate_total_samples(sample_df)
        data_df[_LOCI_COUNT] = self.calculate_total_loci(sample_df)
        data_df[_MUTATION_COUNT] = self.calculate_total_mutations(sample_df)

        style_df = data_df.copy()
        style_df = style_df.applymap(lambda x: "")

        return data_df, style_df

    @staticmethod
    def calculate_total_samples(sample_df):
        sample_df = sample_df.applymap(str)
        sample_df[sample_df == '.'] = None
        sample_df = sample_df.groupby(_GENE_SYMBOL).count()
        sample_df[sample_df > 0] = 1

        return sample_df.apply(sum, 1)

    @staticmethod
    def calculate_total_loci(sample_df):
        return sample_df.groupby(_GENE_SYMBOL).count().ix[:,0]

    @staticmethod
    def calculate_total_mutations(sample_df):
        sample_df = sample_df.applymap(str)
        sample_df[sample_df == '.'] = None
        return sample_df.groupby(_GENE_SYMBOL).count().apply(sum, 1)


class SnpEffFormatRule(object):
    #pylint: disable=too-few-public-methods
    _COLOR_MAP = {"h": "#003366",
                  "m": "#6699FF",
                  "l": "#CCCCFF",
                  "x": "#999999"}

    def format(self, data_df):
        def _determine_format(cell_value):
            for letter in SnpEffFormatRule._COLOR_MAP.keys():
                if letter in cell_value:
                    return letter
            return ""

        data_df = data_df.applymap(str)
        format_df = data_df.applymap(_determine_format)

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
            if len(str(cell_value)) > 0:
                normalized = 100*(cell_value-min_value)/(max_value-min_value)
                return str(int(normalized))
            else:
                return ""

        #TODO: (jebene) make "dbNSFP|overall damaging rank" a constant -pass in?
        format_df =  data_df.drop(self.rank_column, 1)

        min_value = int(min(format_df.min(skipna=True)))
        max_value = int(max(format_df.max(skipna=True)))

        format_df.fillna("", inplace=True)
        format_df = format_df.applymap(_determine_format)

        format_df[self.rank_column] = pd.Series()
        format_df.fillna("", inplace=True)

        for column, dummy in format_df.iteritems():
            for row, dummy in format_df.iterrows():
                styled_cell = self._style(format_df.ix[row, column])
                format_df.ix[row, column] = styled_cell

        return format_df

    @staticmethod
    def _style(cell_value):
        if len(str(cell_value)) > 0:
            cell_value = int(cell_value)
            color = dbNSFPFormatRule._RANGE_RULE[cell_value]
            styled_cell = {"font_size": "12",
                           "bg_color": color.hex,
                           "font_color": dbNSFPFormatRule._FONT_COLOR}
            return styled_cell
        return ""


#TODO: (jebene) hookup rankformatrule
class RankFormatRule(object):
    #pylint: disable=invalid-name, too-few-public-methods
    _RANGE_RULE = list(Color("red").range_to(Color("white"), 101))
    _FONT_COLOR = "#000000"

    def __init__(self):
        #TODO: (jebene) modify regex storage so that annotations and formt rule refer to same format
        self.rank_regex = r".*\|overall .* rank"

    def format(self, data_df):
        def _determine_format(cell_value):
            if len(str(cell_value)) > 0:
                normalized = 100*(cell_value-min_value)/(max_value-min_value)
                return str(int(normalized))
            else:
                return ""

        for col in data_df.columns.values:
            if re.search(self.rank_regex, col):
                min_value = int(data_df[col].min(skipna=True))
                max_value = int(data_df[col].max(skipna=True))

                data_df.fillna("", inplace=True)
                data_df[col] = data_df[col].apply(_determine_format)
            else:
                data_df[col] = data_df[col].apply(lambda x: "")

        for column, dummy in data_df.iteritems():
            for row, dummy in data_df.iterrows():
                styled_cell = self._style(data_df.ix[row, column])
                data_df.ix[row, column] = styled_cell

        return data_df

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

def _combine_dfs(dfs):
    summary_df, dbnsfp_df, snpeff_df = dfs.values()
    combined_df = summary_df.join(dbnsfp_df, how='outer')
    combined_df = combined_df.join(snpeff_df, how='outer')
    combined_df = combined_df.fillna("")
    combined_df.index.names=[_GENE_SYMBOL_OUTPUT_NAME]

    combined_df = combined_df[combined_df.index != "."]

    return combined_df

def _translate_to_excel(data_df, style_df, writer):
#TODO: (jebene): change colors to work in excel
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
    return initial_df.sort(dbNSFP("").damaging_rank_column)

def _rollup(input_file, output_file):
    print "Starting Gene Rollup"

    initial_df = _create_df(input_file)

    annotations = [SummaryColumns(),
                   dbNSFP(dbNSFPFormatRule()),
                   SnpEff(SnpEffFormatRule())]

    annotation_dfs = OrderedDict()
    style_dfs = OrderedDict()

    for annotation in annotations:
        print "Generating {} rollup information".format(annotation.name)
        summarized_df, style_df = annotation.summarize(initial_df)
        annotation_dfs[annotation.name] = summarized_df
        style_dfs[annotation.name] = style_df

    combined_df = _combine_dfs(annotation_dfs)
    combined_df = combined_df.reset_index()
    sorted_df = _sort_by_dbnsfp_rank(combined_df)

    combined_style_df = _combine_dfs(style_dfs)
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
    _SAMPLENAME_REGEX = args.sample_column_regex

    global _GENE_SYMBOL
    _GENE_SYMBOL = args.gene_column_name

    global _XLSX
    _XLSX = args.xlsx

def _add_arg_parse(args):
    parser = argparse.ArgumentParser()
    #pylint: disable=line-too-long
    parser.add_argument("input_file", help=("A tab-delimited file of variants x samples"))
    parser.add_argument("output_file", help=("A tab-delimited file of genes x samples"))
    parser.add_argument("--sample_column_regex", required=True, help="Regex used to define the sample columns in the input file. Default is 'JQ_SUMMARY_SOM_COUNT.*'")
    parser.add_argument("--gene_column_name", required=True, help="Name of gene symbol column in the output file. Default is 'gene symbol'")
    parser.add_argument("--xlsx", action="store_true", help="Write to an xlsx file rather than a tsv file")

    return parser.parse_args(args)

def main():
    args = _add_arg_parse(sys.argv[1:])
    _change_global_variables(args)
    _rollup(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
