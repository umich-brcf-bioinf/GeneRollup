#!/usr/bin/env python
import argparse
import pandas as pd
import re
import sys


_REQUIRED_COLUMNS = set(["GENE_SYMBOL", "dbNSFP_rollup_damaging"])
_SAMPLENAME_REGEX = "JQ_CONS_SOM.*"

def _create_df(input_file):
    initial_df = pd.read_csv(input_file, sep='\t', header=False, dtype='str')
    _validate_df(initial_df)

    return initial_df

def _validate_df(initial_df):
    msg = ("Input file is missing required headers ({})"
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
        return pd.melt(initial_df,
                       id_vars=list(_REQUIRED_COLUMNS),
                       var_name="Sample",
                       value_name="Sample_Data")
    except Exception as excep :
        raise BaseException("Cannot melt dataframe. {0}".format(excep))

def _pivot_df(initial_df):
    initial_df["dbNSFP_rollup_damaging"] = initial_df["dbNSFP_rollup_damaging"].apply(lambda x: int(x))
    return pd.pivot_table(initial_df,
                          index=["GENE_SYMBOL"],
                          columns=["Sample"],
                          values=["dbNSFP_rollup_damaging"],
                          aggfunc=sum)

def _add_arg_parse(args):
    parser = argparse.ArgumentParser()
    #pylint: disable=line-too-long
    parser.add_argument("input_file", help=("A tab-delimited file of variants x samples"))
    parser.add_argument("output_file", help=("A tab-delimited file of genes x samples"))

    return parser.parse_args(args)

def main():
    args = _add_arg_parse(sys.argv[1:])
    initial_df = _create_df(args.input_file)
    condensed_df = _remove_unnecessary_columns(initial_df)
    melted_df = _melt_df(condensed_df)
    pivoted_df = _pivot_df(melted_df)

if __name__ == "__main__":
    main()
