#!/usr/bin/env python
import argparse
import pandas as pd
import re
import sys


_REQUIRED_COLUMNS = set(["GENE_SYMBOL", "dbNSFP_rollup_damaging"])
_SAMPLENAME_REGEX = "JQ_CONS_SOM.*"
_DAMAGING_LABEL = "damaging votes"

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

def _melt_df(initial_df):
    try:
        melted_df = pd.melt(initial_df,
                           id_vars=list(_REQUIRED_COLUMNS),
                           var_name="Sample",
                           value_name="Sample_Data")
        #pylint: disable=line-too-long, unnecessary-lambda
        melted_df["dbNSFP_rollup_damaging"] = melted_df["dbNSFP_rollup_damaging"].apply(lambda x: str(x))
        melted_df["Sample_Data"] = melted_df["Sample_Data"].apply(lambda x: str(x))

        melted_df = melted_df[melted_df["GENE_SYMBOL"] != "."]

        return melted_df
    except Exception as excep :
        raise BaseException("Cannot melt dataframe. {0}".format(excep))

def _pivot_df(initial_df):
    #pylint: disable=line-too-long, unnecessary-lambda
    initial_df["dbNSFP_rollup_damaging"] = initial_df["dbNSFP_rollup_damaging"].apply(lambda x: int(x))
    pivoted_df = pd.pivot_table(initial_df,
                          index=["GENE_SYMBOL"],
                          columns=["Sample"],
                          values=["dbNSFP_rollup_damaging"],
                          aggfunc=sum)
    pivoted_df = pivoted_df.applymap(lambda x: None if x == 0 else str(x))

    return pivoted_df

def _rearrange_columns(initial_df):
    new_column_names = []
    for column_name in list(initial_df.columns.values):
        if type(column_name) is tuple:
            full_source_name, full_sample_name = column_name
            source = full_source_name.split("_")[0]

            patient_prefix = full_sample_name.split("|")[1]
            patient_suffix = full_sample_name.split("|")[2]

            new_column_name = "|".join([source,
                                        _DAMAGING_LABEL,
                                        patient_prefix,
                                        patient_suffix])

        new_column_names.append(new_column_name)
    initial_df.columns = new_column_names

    return initial_df

def _add_arg_parse(args):
    parser = argparse.ArgumentParser()
    #pylint: disable=line-too-long
    parser.add_argument("input_file", help=("A tab-delimited file of variants x samples"))
    parser.add_argument("output_file", help=("A tab-delimited file of genes x samples"))

    return parser.parse_args(args)

def rollup(input_file, output_file):
    initial_df = _create_df(input_file)
    condensed_df = _remove_unnecessary_columns(initial_df)
    melted_df = _melt_df(condensed_df)
    pivoted_df = _pivot_df(melted_df)
    rearranged_df = _rearrange_columns(pivoted_df)

    rearranged_df.to_csv(output_file, sep="\t", index=True)

def main():
    args = _add_arg_parse(sys.argv[1:])

    rollup(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
