#!/usr/bin/env python
import argparse
from collections import OrderedDict
import re
import sys

import pandas as pd


_REQUIRED_COLUMNS = set(["GENE_SYMBOL",
                         "dbNSFP_rollup_damaging",
                         "SNPEFF_TOP_EFFECT_IMPACT"])
_SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
_LOCI_COUNT = "distinct loci"
_SAMPLE_COUNT = "total impacted samples"
_MUTATION_COUNT = "total mutations"
_GENE_SYMBOL = "gene symbol"

class dbNSFP(object):
    #pylint: disable=invalid-name
    def __init__(self):
        self.name = "dbNSFP"
        self.damaging_column = "dbNSFP_rollup_damaging"
        self.damaging_rank_column = "dbNSFP|overall damaging rank"
        self.column_label = "damaging votes"

    def summarize(self, initial_df):
        condensed_df = self._remove_unnecessary_columns(initial_df)
        ranked_df = self._calculate_rank(condensed_df)
        return ranked_df

    def _remove_unnecessary_columns(self, initial_df):
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX)
        required_columns = list(sample_cols.columns)
        required_columns.extend(["GENE_SYMBOL", self.damaging_column])

        condensed_df = pd.DataFrame()
        for col in initial_df.columns:
            if col in required_columns:
                condensed_df[col] = initial_df[col]

        return condensed_df

    def _calculate_rank(self, initial_df):
        #pylint: disable=line-too-long
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values

        initial_df = initial_df.applymap(str)
        initial_df[initial_df == "."] = None

        for sample in sample_cols:
            dummy, samp_number, suffix = sample.split("|")
            sample_name = "|".join([self.name,
                                    self.column_label,
                                    samp_number,
                                    suffix])

            #set sample columns equal to damaging column value
            initial_df[sample][initial_df[sample] != None] = initial_df[self.damaging_column]
            initial_df[sample].fillna(0, inplace=True)
            initial_df[sample_name] = initial_df[sample].apply(int)

            del initial_df[sample]

        ranked_df = initial_df.groupby("GENE_SYMBOL").sum()
        ranked_df = ranked_df.applymap(int)
        ranked_df[self.damaging_rank_column] = ranked_df.apply(sum, 1)

        ranked_df = ranked_df.sort(self.damaging_rank_column, ascending=0)
        ranked_df[self.damaging_rank_column] = ranked_df[self.damaging_rank_column].rank(ascending=0, method="min")

        return ranked_df

    def _change_col_order(self, new_column_names):
        #TODO: (jebene) hookup this method
        ordered_names = []
        samples = ""
        variants = ""
        mutations = ""
        for col in new_column_names:
            if col == _LOCI_COUNT:
                variants = col
            elif col == _SAMPLE_COUNT:
                samples = col
            elif col == _MUTATION_COUNT:
                mutations = col
            elif col == self.damaging_rank_column:
                ordered_names.insert(0, col)
            else:
                ordered_names.append(col)

        if variants:
            ordered_names.insert(0, variants)
        if samples:
            ordered_names.insert(1, samples)
        if mutations:
            ordered_names.insert(2, mutations)
        return ordered_names


class SnpEff(object):
    _RANK_SCORES =  {"HIGH": 100000.0,
                     "MODERATE": 1,
                     "LOW": 1/100000.0,
                     "MODIFIER": 10**12}
    _RANK_ABBREVS =  {"HIGH": "h", "MODERATE": "m", "LOW": "l", "MODIFIER": "x"}

    #pylint: disable=invalid-name
    def __init__(self):
        self.name = "SnpEff"
        self.impact_column = "SNPEFF_TOP_EFFECT_IMPACT"
        self.impact_rank_column = "SnpEff|overall impact rank"
        self.impact_score = "impact score"
        self.impact_score_column = "SnpEff|overall impact score"
        self.column_label = "impact"
        self.impact_category = "SnpEff|impact category|{}"

    def summarize(self, initial_df):
        condensed_df = self._remove_unnecessary_columns(initial_df)
#         modified_df, score = self._preprocess_df(condensed_df)
        ranked_df = self._calculate_rank(condensed_df)
        return ranked_df

    def _remove_unnecessary_columns(self, initial_df):
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX)
        required_columns = list(sample_cols.columns)
        required_columns.extend(["GENE_SYMBOL", self.impact_column])

        condensed_df = pd.DataFrame()
        for col in initial_df.columns:
            if col in required_columns:
                condensed_df[col] = initial_df[col]

        return condensed_df

    def _calculate_rank(self, initial_df):
        #TODO: (jebene) add counts for HIGH, MODERATE, etc.

        #pylint: disable=line-too-long
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values

        initial_df = initial_df.applymap(str)
        initial_df[initial_df == "."] = None

        scored_df = pd.DataFrame()
        scored_df["GENE_SYMBOL"] = initial_df["GENE_SYMBOL"]

        for sample in sample_cols:
            dummy, samp_number, suffix = sample.split("|")
            sample_name = "|".join([self.name,
                                    self.column_label,
                                    samp_number,
                                    suffix])

            #set sample columns equal to impact column value
            initial_df[sample][initial_df[sample] != None] = initial_df[self.impact_column]
            initial_df[sample_name] = initial_df[sample].map(SnpEff._RANK_ABBREVS)
            scored_df[sample + "_SCORE"] = initial_df[sample].map(SnpEff._RANK_SCORES)

            del initial_df[sample]

        score = scored_df.groupby("GENE_SYMBOL").sum().apply(sum, 1)
        ranked_df = initial_df.groupby("GENE_SYMBOL").sum()

        if self.impact_column in ranked_df.columns.values:
            del ranked_df[self.impact_column]

#         ranked_df["count"] = ranked_df.apply(sum, 1)
#         print ranked_df

        ranked_df[self.impact_score_column] = score

        ranked_df = ranked_df.sort(self.impact_score_column, ascending=0)
        ranked_df[self.impact_rank_column] = ranked_df[self.impact_score_column].rank(ascending=0, method="min")

        return ranked_df

    def _change_col_order(self, new_column_names):
        #TODO: (jebene) hookup this method
        rank = []
        score = []
        category = []
        impact = []
        for col in new_column_names:
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

        return ordered_names


class SummaryColumns(object):
    def __init__(self):
        self.name = "Summary Columns"

    def summarize(self, initial_df):
        sample_df =  initial_df.filter(regex=_SAMPLENAME_REGEX)
        sample_df["GENE_SYMBOL"] = initial_df["GENE_SYMBOL"]

        summary_df = pd.DataFrame()
        summary_df[_SAMPLE_COUNT] = self.calculate_total_samples(sample_df)
        summary_df[_LOCI_COUNT] = self.calculate_total_loci(sample_df)
        summary_df[_MUTATION_COUNT] = self.calculate_total_mutations(sample_df)

        return summary_df

    @staticmethod
    def calculate_total_samples(sample_df):
        sample_df = sample_df.applymap(str)
        sample_df[sample_df == '.'] = None
        sample_df = sample_df.groupby("GENE_SYMBOL").count()
        sample_df[sample_df > 0] = 1

        return sample_df.apply(sum, 1)

    @staticmethod
    def calculate_total_loci(sample_df):
        return sample_df.groupby("GENE_SYMBOL").count().ix[:,0]

    @staticmethod
    def calculate_total_mutations(sample_df):
        sample_df = sample_df.applymap(str)
        sample_df[sample_df == '.'] = None
        return sample_df.groupby("GENE_SYMBOL").count().apply(sum, 1)

    @staticmethod
    def rearrange_columns(initial_df):
        initial_df.index.names=[_GENE_SYMBOL]
        return initial_df

def _create_df(input_file):
    initial_df = pd.read_csv(input_file, sep='\t', header=False, dtype='str')
    _validate_df(initial_df)

    return initial_df

def _validate_df(initial_df):
    header = set(initial_df.columns.values)
    missing_columns = _REQUIRED_COLUMNS.difference(header)
    msg = ("Input file is missing required headers ({}). "
           "Review input and try again."
           ).format(missing_columns)
    if missing_columns:
        raise BaseException(msg)

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
    combined_df = combined_df.fillna(0)

    combined_df = combined_df.fillna("")

    return combined_df

def _rollup(input_file, output_file):
    print "Starting Gene Rollup"

    initial_df = _create_df(input_file)

    annotations = [SummaryColumns(), dbNSFP(), SnpEff()]
    annotation_dfs = OrderedDict()

    for annotation in annotations:
        print "Generating {} rollup information".format(annotation.name)
        summarized_df = annotation.summarize(initial_df)
        annotation_dfs[annotation.name] = summarized_df

    combined_df = _combine_dfs(annotation_dfs)
    combined_df.to_csv(output_file, sep="\t", index=True)

    print "Done."

def _add_arg_parse(args):
    parser = argparse.ArgumentParser()
    #pylint: disable=line-too-long
    parser.add_argument("input_file", help=("A tab-delimited file of variants x samples"))
    parser.add_argument("output_file", help=("A tab-delimited file of genes x samples"))

    return parser.parse_args(args)

def main():
    args = _add_arg_parse(sys.argv[1:])
    _rollup(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
