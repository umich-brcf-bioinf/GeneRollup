#!/usr/bin/env python
import argparse
from collections import OrderedDict
import re
import sys

import numpy as np
import pandas as pd


_REQUIRED_COLUMNS = set(["GENE_SYMBOL", #chrom pos ref alt
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
        melted_df = self.melt_df(initial_df)
        pivoted_df = self.pivot_df(melted_df)
        ranked_df = self.calculate_rank(pivoted_df)
        return ranked_df

    def melt_df(self, initial_df):
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
        try:
            melted_df = pd.melt(initial_df,
                               id_vars=list(_REQUIRED_COLUMNS),
                               value_vars=list(sample_cols),
                               var_name="Sample",
                               value_name="Sample_Data")

            #pylint: disable=line-too-long, unnecessary-lambda
            melted_df[self.damaging_column] = melted_df[self.damaging_column].apply(lambda x: str(x))
            melted_df = melted_df.applymap(lambda x: str(x))
            melted_df = melted_df[melted_df["GENE_SYMBOL"] != "."]

            return melted_df

        except Exception as excep :
            raise BaseException("Cannot melt dataframe. {0}".format(excep))

    def pivot_df(self, initial_df):
        #pylint: disable=line-too-long, unnecessary-lambda
        initial_df[self.damaging_column] = initial_df[self.damaging_column].apply(lambda x: int(x))
        pivoted_df = pd.pivot_table(initial_df,
                              index=["GENE_SYMBOL"],
                              columns=["Sample"],
                              values=[self.damaging_column],
                              aggfunc=sum)
        pivoted_df = pivoted_df.applymap(lambda x: None if x == 0 else str(x))
        pivoted_df.fillna(value="0", inplace=True)

        return pivoted_df

    def calculate_rank(self, initial_df):
        #pylint: disable=line-too-long, unnecessary-lambda
        initial_df[self.damaging_rank_column] = initial_df[self.damaging_column].sum(axis=1)
        initial_df = initial_df.sort(self.damaging_rank_column, ascending=0)
        initial_df[self.damaging_rank_column] = initial_df[self.damaging_rank_column].rank(ascending=0, method="min")

        try:
            initial_df[self.damaging_rank_column] = initial_df[self.damaging_rank_column].apply(lambda x: str(int(x)))
        except ValueError:
            pass

        return initial_df

    def rearrange_columns(self, initial_df):
        new_column_names = []
        totals = []

        for column_name in list(initial_df.columns.values):
            if type(column_name) is tuple:
                if len(column_name[1]) > 0:
                    if column_name == (self.damaging_column, _SAMPLE_COUNT):
                        totals.append(column_name[1])
                    else:
                        full_source_name, full_sample_name = column_name
                        source = full_source_name.split("_")[0]

                        patient_prefix = full_sample_name.split("|")[1]
                        patient_suffix = full_sample_name.split("|")[2]

                        new_column_name = "|".join([source,
                                                    self.column_label,
                                                    patient_prefix,
                                                    patient_suffix])
                        new_column_names.append(new_column_name)
                else:
                    new_column_names.append(column_name[0])

        all_column_names = totals + new_column_names

        initial_df.columns = all_column_names
        ordered_names = self._change_col_order(all_column_names)

        initial_df = initial_df[ordered_names]
        initial_df.index.names=[_GENE_SYMBOL]

        return initial_df

    def _change_col_order(self, new_column_names):
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
    _POSSIBLE_VALUES = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
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
        melted_df = self.melt_df(initial_df)
        pivoted_df = self.pivot_df(melted_df)
        ranked_df = self.calculate_rank(pivoted_df)
        return ranked_df

    def melt_df(self, initial_df):
        sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX).columns.values
        try:
            melted_df = pd.melt(initial_df,
                               id_vars=list(_REQUIRED_COLUMNS),
                               value_vars=list(sample_cols),
                               var_name="Sample",
                               value_name="Sample_Data")
            #pylint: disable=line-too-long, unnecessary-lambda
            melted_df[self.impact_column] = melted_df[self.impact_column].apply(lambda x: str(x))
            melted_df = melted_df.applymap(lambda x: str(x))

            melted_df = melted_df[melted_df["GENE_SYMBOL"] != "."]

            return melted_df

        except Exception as excep :
            raise BaseException("Cannot melt dataframe. {0}".format(excep))

    def pivot_df(self, initial_df):
        #pylint: disable=line-too-long, unnecessary-lambda
        pivoted_df = pd.pivot_table(initial_df,
                              index=["GENE_SYMBOL", "Sample"],
                              columns=[self.impact_column],
                              values=["Sample_Data"],
                              #pylint: disable=no-member
                              aggfunc=np.count_nonzero,
                              fill_value=0)
        pivoted_df.fillna(value="0", inplace=True)
        pivoted_df = pivoted_df.applymap(lambda x: int(x))

        pivoted_df = pivoted_df["Sample_Data"]

        return pivoted_df

    def calculate_rank(self, initial_df):
        #pylint: disable=line-too-long
        for item in SnpEff._POSSIBLE_VALUES:
            if item not in initial_df.columns:
                initial_df[item] = 0
            initial_df[item + "_initial_sum"] = initial_df[item].map(int)

        impact_df = self._calculate_impacts(initial_df)
        scored_df = self._calculate_score(impact_df)

        expanded_df = scored_df.unstack()

        expanded_df[self.impact_score_column] = expanded_df[self.impact_score].sum(axis=1)
        expanded_df = self._calculate_impact_score(expanded_df)
        expanded_df = self._calculate_impact_rank(expanded_df)

        del expanded_df[self.impact_score]
        for item in SnpEff._POSSIBLE_VALUES:
            del expanded_df[item + "_initial_sum"]

        try:
            expanded_df[self.impact_rank_column] = expanded_df[self.impact_rank_column].apply(lambda x: str(int(x)))
        except ValueError:
            pass

        return expanded_df

    def _calculate_impacts(self, initial_df):
        high = initial_df["HIGH"].apply(lambda x: "h" * x)
        moderate =  initial_df["MODERATE"].apply(lambda x: "m" * x)
        low = initial_df["LOW"].apply(lambda x: "l" * x)
        modifier = initial_df["MODIFIER"].apply(lambda x: "x" * x)

        initial_df[self.name] = high + moderate + low + modifier

        return initial_df

    def _calculate_score(self, initial_df):
        high = initial_df["HIGH"] * 100000.0
        moderate = initial_df["MODERATE"]
        low = initial_df["LOW"]/100000.0
        modifier = initial_df["MODIFIER"]/10**12

        initial_df[self.impact_score] = high + moderate + low + modifier

        del initial_df["HIGH"]
        del initial_df["MODERATE"]
        del initial_df["LOW"]
        del initial_df["MODIFIER"]

        return initial_df

    def _calculate_impact_score(self, initial_df):
        #pylint: disable=line-too-long
        for item in SnpEff._POSSIBLE_VALUES:
            initial_df[self.impact_category.format(item)] = initial_df[item + "_initial_sum"].sum(axis=1)

        return initial_df

    def _calculate_impact_rank(self, initial_df):
        #pylint: disable=line-too-long
        initial_df = initial_df.sort(self.impact_score_column, ascending=0)
        initial_df[self.impact_rank_column] = initial_df[self.impact_score_column].rank(ascending=0, method="min")

        return initial_df

    def rearrange_columns(self, initial_df):
        new_column_names = []
        totals = []

        for column_name in list(initial_df.columns.values):
            if type(column_name) is tuple:
                if len(column_name[1]) > 0:
                    if column_name == (self.name , _SAMPLE_COUNT):
                        del initial_df[column_name]
                    else:
                        full_source_name, full_sample_name = column_name
                        source = full_source_name.split("_")[0]

                        patient_prefix = full_sample_name.split("|")[1]
                        patient_suffix = full_sample_name.split("|")[2]

                        new_column_name = "|".join([source,
                                                    self.column_label,
                                                    patient_prefix,
                                                    patient_suffix])
                        new_column_names.append(new_column_name)
                else:
                    new_column_names.append(column_name[0])

        all_column_names = totals + new_column_names

        initial_df.columns = all_column_names
        ordered_names = self._change_col_order(all_column_names)

        initial_df = initial_df[ordered_names]
        initial_df.index.names=[_GENE_SYMBOL]

        return initial_df

    def _change_col_order(self, new_column_names):
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
        self.groupby_column = "GENE_SYMBOL"

    def summarize(self, initial_df):
        sample_df =  initial_df.filter(regex=_SAMPLENAME_REGEX)
        sample_df[self.groupby_column] = initial_df[self.groupby_column]

        summary_df = pd.DataFrame()
        summary_df[_SAMPLE_COUNT] = self.calculate_total_samples(sample_df)
        summary_df[_LOCI_COUNT] = self.calculate_total_loci(sample_df)
        summary_df[_MUTATION_COUNT] = self.calculate_total_mutations(sample_df)

        return summary_df

    def calculate_total_samples(self, sample_df):
        sample_df = sample_df.applymap(str)
        sample_df[sample_df == '.'] = None
        sample_df = sample_df.groupby(self.groupby_column).count()
        sample_df[sample_df > 0] = 1

        return sample_df.apply(sum, 1)

    def calculate_total_loci(self, sample_df):
        return sample_df.groupby(self.groupby_column).count().ix[:,0]

    def calculate_total_mutations(self, sample_df):
        sample_df = sample_df.applymap(str)
        sample_df[sample_df == '.'] = None
        return sample_df.groupby(self.groupby_column).count().apply(sum, 1)

    @staticmethod
    def rearrange_columns(initial_df):
        initial_df.index.names=[_GENE_SYMBOL]
        return initial_df

def _create_df(input_file):
    initial_df = pd.read_csv(input_file, sep='\t', header=False, dtype='str')
    _validate_df(initial_df)

    return initial_df

def _validate_df(initial_df):
    msg = ("Input file is missing required headers ({}). "
           "Review input and try again.",
           _REQUIRED_COLUMNS)
    header = set(initial_df.columns.values)
    if not _REQUIRED_COLUMNS.issubset(header):
        raise BaseException(msg)

    sample_column = 0
    for column in header:
        if re.search(_SAMPLENAME_REGEX, column):
            sample_column = 1
            break
    if not sample_column:
        raise BaseException(msg)

def _remove_unnecessary_columns(initial_df):
    sample_cols =  initial_df.filter(regex=_SAMPLENAME_REGEX)
    required_columns = list(sample_cols.columns)
    required_columns.extend(_REQUIRED_COLUMNS)

    for col in initial_df.columns:
        if col not in required_columns:
            del initial_df[col]

    return initial_df

def _combine_dfs(dfs):
    summary_df, dbnsfp_df, snpeff_df = dfs.values()
    combined_df = summary_df.join(dbnsfp_df, how='outer')
    combined_df = combined_df.join(snpeff_df, how='outer')
    combined_df = combined_df.fillna(0)

    combined_df = combined_df.fillna("")

    return combined_df

def _add_arg_parse(args):
    parser = argparse.ArgumentParser()
    #pylint: disable=line-too-long
    parser.add_argument("input_file", help=("A tab-delimited file of variants x samples"))
    parser.add_argument("output_file", help=("A tab-delimited file of genes x samples"))

    return parser.parse_args(args)

def _rollup(input_file, output_file):
    print "Starting Gene Rollup"

    initial_df = _create_df(input_file)
    condensed_df = _remove_unnecessary_columns(initial_df)

    annotations = [SummaryColumns(), dbNSFP(), SnpEff()]
    annotation_dfs = OrderedDict()

    for annotation in annotations:
        print "Generating {} rollup information".format(annotation.name)
        summarized_df = annotation.summarize(condensed_df)
        rearranged_df = annotation.rearrange_columns(summarized_df)
        annotation_dfs[annotation.name] = rearranged_df

    combined_df = _combine_dfs(annotation_dfs)
    combined_df.to_csv(output_file, sep="\t", index=True)

    print "Done."


def main():
    args = _add_arg_parse(sys.argv[1:])
    _rollup(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
