#!/usr/bin/env python
from StringIO import StringIO
import argparse
import pandas as pd
import re
import sys


_REQUIRED_COLUMNS = set(["GENE_SYMBOL", "dbNSFP_rollup_damaging"])
_SAMPLENAME_REGEX = "JQ_CONS_SOM.*"

def _create_df(input_file):
    initial_df = pd.read_csv(input_file, sep='\t', header=False, dtype='str')
    
    return initial_df

def _validate_df(df):
    msg = ("Input file is missing required headers ({}). Review input and try again.",
           _REQUIRED_COLUMNS)
    header = set(df.columns.values)
    if not _REQUIRED_COLUMNS.issubset(header):
        raise BaseException(msg)

    sample_column = 0
    for column in header:
        if re.search(_SAMPLENAME_REGEX, column):
            sample_column = 1
            break
    if not sample_column:
        raise BaseException(msg)

def _remove_unnecessary_columns(df):
    sample_cols =  df.filter(regex=_SAMPLENAME_REGEX)
    required_columns = list(sample_cols.columns)
    required_columns.extend(_REQUIRED_COLUMNS)

    for col in df.columns:
        if col not in required_columns:
            del df[col]

    return df

def _add_arg_parse(args):
    parser = argparse.ArgumentParser()
    #pylint: disable=line-too-long
    parser.add_argument("input_file", help=("A tab-delimited file of variants x samples"))
    parser.add_argument("output_file", help=("A tab-delimited file of genes x samples"))

    return parser.parse_args(args)

def main():
    args = _add_arg_parse(sys.argv[1:])
    initial_df = _create_df(args.input_file)
    _validate_df(initial_df)

if __name__ == "__main__":
    main()