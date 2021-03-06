#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=maybe-no-member, too-few-public-methods, no-member
from __future__ import absolute_import

from argparse import Namespace
from collections import OrderedDict
import filecmp
from io import StringIO
import os
import unittest

import numpy
import pandas as pd

from testfixtures import TempDirectory

import generollup.rollup as rollup

#TODO: cgates: How about we fix the chained assigments and get rid of this warning suppression.
pd.set_option('mode.chained_assignment', None)

def dataframe(input_data, sep="|", dtype=None):
    return pd.read_csv(StringIO(input_data), sep=sep, header=0, dtype=dtype)

class MockFormatRule(object):
    def __init__(self, format_df):
        self.format_df = format_df
        self.last_data_df = None
        self.last_format_df = None

    def format(self, data_df):
        self.last_data_df = data_df
        return self.format_df

    def style(self, format_df):
        self.last_format_df = format_df
        return self.format_df

class GeneRollupTestCase(unittest.TestCase):
    def setUp(self):
        rollup._DBNSFP_COLUMN = 'dbNSFP_rollup_damaging'
        rollup._EFFECT_COLUMN = 'SNPEFF_TOP_EFFECT_IMPACT'
        rollup._GENE_SYMBOL = 'GENE_SYMBOL'

    def test_create_df(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|sampleA
1\t2\t3\t4\t5\t6\t7\t8'''

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        actual_df = rollup._create_df(StringIO(input_string), args)
        self.assertEquals(["GENE_SYMBOL", "dbNSFP_rollup_damaging", "SNPEFF_TOP_EFFECT_IMPACT", "JQ_SUMMARY_SOM_COUNT|sampleA"],
                            list(actual_df.columns.values))

    def test_create_df_missingDbnsfpAndSnpeffOkay(self):
        input_string =\
'''GENE_SYMBOL\theaderA\theaderB\tJQ_SUMMARY_SOM_COUNT|sampleA
foo\t1\t2\t0'''

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        actual_df = rollup._create_df(StringIO(input_string), args)
        self.assertEquals(["GENE_SYMBOL", "headerA", "headerB", "JQ_SUMMARY_SOM_COUNT|sampleA"],
                            list(actual_df.columns.values))

    def test_create_df_missingDbNsfpOkay(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|sampleA|TUMOR
1\t2\t3'''

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name='SNPEFF_TOP_EFFECT_IMPACT',
                         tsv=False)

        actual_df = rollup._create_df(StringIO(input_string), args)
        self.assertEquals(["GENE_SYMBOL", "SNPEFF_TOP_EFFECT_IMPACT", "JQ_SUMMARY_SOM_COUNT|sampleA|TUMOR"],
                            list(actual_df.columns.values))

    def test_create_df_missingEffectOkay(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|sampleA|TUMOR
1\t2\t3'''

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         dbnsfp_column_name='dbNSFP_rollup_damaging',
                         tsv=False)

        actual_df = rollup._create_df(StringIO(input_string), args)
        self.assertEquals(["GENE_SYMBOL", "dbNSFP_rollup_damaging", "JQ_SUMMARY_SOM_COUNT|sampleA|TUMOR"],
                            list(actual_df.columns.values))

    def test_create_df_missingSamples(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tSNPEFF_TOP_EFFECT_IMPACT
1\t2\t3
1\t2\t3'''

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(P.)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         dbnsfp_column_name='dbNSFP_rollup_damaging',
                         effect_column_name='SNPEFF_TOP_EFFECT_IMPACT',
                         tsv=False)

        self.assertRaisesRegexp(rollup.UsageError,
                                "Cannot determine sample genotype columns with supplied regex.*",
                                rollup._create_df,
                                StringIO(input_string),
                                args)

    #TODO: (jebene) I can't figure out how to initialize this as having null values
    def xtest_create_df_removesIntergenicVariants(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT
BRCA1\t2\t3
0\t4\t5
6\t7\t8'''
        input_string = input_string.replace("0", numpy.nan)
        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name="GENE_SYMBOL",
                         tsv=False)

        actual_df = rollup._create_df(StringIO(input_string), args)

        self.assertEquals(["GENE_SYMBOL", "SNPEFF_TOP_EFFECT_IMPACT", "JQ_SUMMARY_SOM_COUNT"],
                            list(actual_df.columns.values))
        self.assertEquals(["BRCA1", "6"], list(actual_df["GENE_SYMBOL"].values))

    def test_sort_by_dbnsfp_rank(self):
        input_string =\
'''gene_symbol\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\teffect_annotation|overall_effect_rank\tdbNSFP_annotation|overall_damaging_rank
BRCA1\th\t7\t2
EGFR\tm\t4\t3
SON\tm\t5\t1
BRCA2\tm\t5\t1
CREBBP\thhh\t6\t1'''
        input_df = dataframe(input_string, sep="\t")
        sorted_df = rollup._sort_by_dbnsfp_rank(input_df)

        self.assertEquals(["BRCA2", "SON", "CREBBP", "BRCA1", "EGFR"], list(sorted_df["gene_symbol"].values))
        self.assertEquals([1, 1, 1, 2, 3], list(sorted_df["dbNSFP_annotation|overall_damaging_rank"].values))

    def test_combine_dfs(self):
        summary_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL
BRCA1\t1
EGFR\t1
CREBBP\t1'''
        summary_df = dataframe(summary_string, sep="\t")
        summary_df = summary_df.set_index(["GENE_SYMBOL"])

        dbNSFP_string =\
'''GENE_SYMBOL\tdbNSFP|P1|NORMAL
BRCA1\t2
CREBBP\t4'''
        dbNSFP_df = dataframe(dbNSFP_string, sep="\t")
        dbNSFP_df = dbNSFP_df.set_index(["GENE_SYMBOL"])

        snpEff_string =\
'''GENE_SYMBOL\tsnpEff|P1|NORMAL
BRCA1\th
CREBBP\thh'''
        snpEff_df = dataframe(snpEff_string, sep="\t")
        snpEff_df = snpEff_df.set_index(["GENE_SYMBOL"])

        dfs = OrderedDict()
        dfs["summary"] = summary_df
        dfs["dbNSFP"] = dbNSFP_df
        dfs["snpEff"] = snpEff_df

        actual = rollup._combine_dfs(dfs)

        self.assertEquals(["BRCA1", "CREBBP"], list(actual.index.values))

    def xtest_translate_to_excel(self):
        with TempDirectory() as output_dir:
            output_dir.write("output.xlsx", "")
            output_file = os.path.join(output_dir.path, "output.xlsx")
            data_string =\
    '''gene symbol|PATIENT_A_SnpEff|PATIENT_A_dbNSFP|SnpEff_overall_impact_rank
    MOD|mml|12|2
    NULL1|||
    HIGH|hhmlx|4|1'''
            data_df = dataframe(data_string)
            data_df.fillna("", inplace=True)

            style_string = \
    '''gene symbol|PATIENT_A_SnpEff|PATIENT_A_dbNSFP|SnpEff_overall_impact_rank
    MOD|||
    HIGH|||
    NULL1|||'''

            style_df = dataframe(style_string)
            style_df["PATIENT_A_SnpEff"] = [{"font_size": "4", "bg_color": "#6699FF", "font_color": "#6699FF"},
                                            "",
                                            {"font_size": "4", "bg_color": "#003366", "font_color": "#003366"}]
            style_df["PATIENT_A_dbNSFP"] = [{"font_size": "12", "bg_color": "#ffa500", "font_color": "#000000"},
                                            "",
                                            {"font_size": "12", "bg_color": "white", "font_color": "#003366"}]
            style_df["SnpEff_overall_impact_rank"] = [{"font_size": "12", "bg_color": "white", "font_color": "#000000"},
                                                      "",
                                                      {"font_size": "12", "bg_color": "red", "font_color": "#000000"}]
            style_df.fillna("", inplace=True)

            writer = pd.ExcelWriter(output_file, engine="xlsxwriter")
            rollup._translate_to_excel(data_df, style_df, writer)

            script_dir = os.path.dirname(os.path.realpath(__file__))
            expected_output = os.path.join(script_dir,
                                           "functional_tests",
                                           "translate_to_excel",
                                           "expected_output.xlsx")
            self.assertEquals(True, filecmp.cmp(expected_output, output_file))

    def test_reset_style_gene_values(self):
        data_string =\
'''gene_symbol|PATIENT_A_SnpEff
BRCA1|{"foo": "bar"}
TANK|{"foo": "bar"}
CREBBP|{"foo": "bar"}'''
        data_df = dataframe(data_string)
        actual = rollup._reset_style_gene_values(data_df)
        actual = actual.applymap(str)

        expected_string =\
'''gene_symbol|PATIENT_A_SnpEff
{}|{"foo": "bar"}
{}|{"foo": "bar"}
{}|{"foo": "bar"}'''
        expected = dataframe(expected_string)

        self.assertEquals(list(expected["gene_symbol"].values),
                          list(actual["gene_symbol"].values))

        self.assertEquals(list(expected["PATIENT_A_SnpEff"].values),
                          list(actual["PATIENT_A_SnpEff"].values))


class dbNSFPTestCase(unittest.TestCase):
    def setUp(self):
        rollup._SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
        rollup._GENE_SYMBOL = "GENE_SYMBOL"
        rollup._XLSX = False

    def tearDown(self):
        pass

    def test_remove_unnecessary_columns(self):
        FORMAT_DF = pd.DataFrame([[42] * 4] * 1)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name="GENE_SYMBOL",
                         tsv=False)

        dbNSFP = rollup.dbNSFP([formatRule], args)

        input_string = \
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tBAZ\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR\tFOO\tBAR'''
        input_df = dataframe(input_string, sep="\t")

        actual = dbNSFP._remove_unnecessary_columns(input_df)

        self.assertEquals(4, len(actual.columns))
        self.assertNotIn("BAZ", actual.columns)
        self.assertNotIn("FOO", actual.columns)
        self.assertNotIn("BAR", actual.columns)

    def test_remove_invalid_rows(self):
        FORMAT_DF = pd.DataFrame([[42] * 4] * 2)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         )

        dbNSFP = rollup.dbNSFP([formatRule], args)

        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\t2\t1\t1
BRCA1\t0\t.\t1
BRCA1\t\t1\t1
BRCA1\t.\t1\t1
CREBBP\t3\t0\t.'''

        input_df = dataframe(input_string, sep="\t")
        actual = dbNSFP._remove_invalid_rows(input_df)

        self.assertEquals(["2", "3"], list(actual["dbNSFP_rollup_damaging"].values))
        self.assertEquals(["1", "0"], list(actual["JQ_SUMMARY_SOM_COUNT|P1|TUMOR"].values))
        self.assertEquals(["1", "."], list(actual["JQ_SUMMARY_SOM_COUNT|P2|TUMOR"].values))

    def test_summarize_dataMatrix(self):
        FORMAT_DF = pd.DataFrame([[42] * 4] * 2)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name=rollup._GENE_SYMBOL,
                         input_gene_column_name=rollup._GENE_SYMBOL,
                         tsv=False)

        dbNSFP = rollup.dbNSFP([formatRule], args)

        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\t2\t1\t1
BRCA1\t5\t.\t1
CREBBP\t3\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        input_df = input_df.applymap(str)
        (data_df, style_dfs) = dbNSFP.summarize(input_df)

        self.assertEquals(1, len(style_dfs))
        self.assertEquals(data_df.shape, style_dfs[0].shape)

        data_df = data_df.applymap(str)
        expected_string =\
'''GENE_SYMBOL\tdbNSFP_annotation|overall_damaging_rank\tdbNSFP_annotation|damaging_total\tdbNSFP_annotation|damaging_votes|P1\tdbNSFP_annotation|damaging_votes|P2
BRCA1\t1\t9\t2\t7
CREBBP\t2\t0\t\t'''
        expected_df = dataframe(expected_string, sep="\t", dtype=str)
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals('\t'.join(expected_df.columns.values), '\t'.join(data_df.columns.values))
        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in data_df.values])

    def test_summarize_dataMatrixIgnoresNullOrZeroDamagingCounts(self):
        FORMAT_DF = pd.DataFrame([[42] * 4] * 1)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name=rollup._GENE_SYMBOL,
                         input_gene_column_name=rollup._GENE_SYMBOL,
                         tsv=False)

        dbNSFP = rollup.dbNSFP([formatRule], args)

        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\t0\t1\t1
BRCA1\t1\t.\t1
CREBBP\t.\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        input_df = input_df.applymap(str)
        (data_df, style_dfs) = dbNSFP.summarize(input_df)

        self.assertEquals(1, len(style_dfs))
        self.assertEquals(data_df.shape, style_dfs[0].shape)

        data_df = data_df.applymap(str)
        expected_string =\
'''GENE_SYMBOL\tdbNSFP_annotation|overall_damaging_rank\tdbNSFP_annotation|damaging_total\tdbNSFP_annotation|damaging_votes|P1\tdbNSFP_annotation|damaging_votes|P2
BRCA1\t1\t1\t\t1'''
        expected_df = dataframe(expected_string, sep="\t", dtype=str)
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals('\t'.join(expected_df.columns.values), '\t'.join(data_df.columns.values))
        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in data_df.values])

    def test_summarize_formatMatrix(self):
        FORMAT_DF = pd.DataFrame([[42] * 4] * 2)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name=rollup._GENE_SYMBOL,
                         input_gene_column_name=rollup._GENE_SYMBOL,
                         tsv=False)

        dbNSFP = rollup.dbNSFP([formatRule], args)

        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\t2\t1\t1
BRCA1\t5\t.\t1
CREBBP\t3\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        input_df = input_df.applymap(str)
        (dummy, format_dfs) = dbNSFP.summarize(input_df)

        self.assertEquals(1, len(format_dfs))
        self.assertIs(FORMAT_DF, format_dfs[0])

    def test_summarize_multipleFormatRules(self):
        FORMAT_DF1 = pd.DataFrame([[42] * 4] * 2)
        formatRule1 = MockFormatRule(FORMAT_DF1)

        FORMAT_DF2 = pd.DataFrame([["A"] * 4] * 2)
        formatRule2 = MockFormatRule(FORMAT_DF2)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        dbNSFP = rollup.dbNSFP([formatRule1, formatRule2], args)

        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t2\t1\t1
BRCA1\t5\t.\t1
CREBBP\t3\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        input_df = input_df.applymap(str)
        (data_df, format_dfs) = dbNSFP.summarize(input_df)

        self.assertIs(data_df, formatRule1.last_data_df)
        self.assertEquals(2, len(format_dfs))
        self.assertIs(FORMAT_DF1, format_dfs[0])
        self.assertIs(FORMAT_DF2, format_dfs[1])

    def test_build_damaging_votes(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\t2\t1\t1
BRCA1\t3\t.\t1
BRCA1\t3\t\t1
CREBBP\t2\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         dbnsfp_column_name='dbNSFP_rollup_damaging',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        dbNSFP = rollup.dbNSFP([formatRule], args)

        damaging_votes_df = dbNSFP._build_gene_sample_damaging_votes(input_df)
        self.assertEquals(["BRCA1", "CREBBP"], list(damaging_votes_df.index.values))
        import numpy as np
        self.assertEquals([2.0, '.'], list(damaging_votes_df["dbNSFP_annotation|damaging_votes|P1"].values))
        self.assertTrue([8.0, '.'], list(damaging_votes_df["dbNSFP_annotation|damaging_votes|P2"].values))


    def test_build_damaging_ranks(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t2\t1\t1
BRCA1\t3\t.\t1
CREBBP\t2\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         dbnsfp_column_name='dbNSFP_rollup_damaging',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        dbNSFP = rollup.dbNSFP([formatRule], args)

        damaging_votes_df = dbNSFP._build_gene_sample_damaging_votes(input_df)
        ranked_df = dbNSFP._build_damaging_ranks(damaging_votes_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 2], list(ranked_df["dbNSFP_annotation|overall_damaging_rank"].values))

    def test_build_damaging_ranks_tie(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
TANK\t2\t.\t1
BRCA1\t5\t.\t1
EGFR\t5\t.\t1
CREBBP\t5\t1\t.
BRCA2\t13\t.\t1'''
        input_df = dataframe(input_string, sep="\t")
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         dbnsfp_column_name='dbNSFP_rollup_damaging',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        dbNSFP = rollup.dbNSFP([formatRule], args)

        damaging_votes_df = dbNSFP._build_gene_sample_damaging_votes(input_df)
        ranked_df = dbNSFP._build_damaging_ranks(damaging_votes_df)
        self.assertEquals(["BRCA1", "BRCA2", "CREBBP", "EGFR", "TANK"], list(ranked_df.index.values))
        self.assertEquals([2, 1, 2, 2, 5], list(ranked_df["dbNSFP_annotation|overall_damaging_rank"].values))

class EffectTestCase(unittest.TestCase):
    def setUp(self):
        rollup._SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
        rollup._GENE_SYMBOL = "GENE_SYMBOL"
        rollup._XLSX = False

    def tearDown(self):
        pass

    def test_remove_unnecessary_columns(self):
        FORMAT_DF = pd.DataFrame([[42] * 4] * 1)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name = 'SNPEFF_TOP_EFFECT_IMPACT',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        snpEff = rollup.Effect([formatRule], args)

        input_string = \
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tBAZ\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR\tFOO\tBAR'''
        input_df = dataframe(input_string, sep="\t")

        actual = snpEff._remove_unnecessary_columns(input_df)

        self.assertEquals(4, len(actual.columns))
        self.assertNotIn("BAZ", actual.columns)
        self.assertNotIn("FOO", actual.columns)
        self.assertNotIn("BAR", actual.columns)

    def test_remove_invalid_rows(self):
        FORMAT_DF = pd.DataFrame([[42] * 4] * 2)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name = 'SNPEFF_TOP_EFFECT_IMPACT',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        effect = rollup.Effect([formatRule], args)

        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\tHIGH\t1\t1
BRCA1\tLOW\t.\t1
BRCA1\t\t1\t1
BRCA1\t.\t1\t1
CREBBP\tMODIFIER\t0\t.'''

        input_df = dataframe(input_string, sep="\t")
        actual = effect._remove_invalid_rows(input_df)

        self.assertEquals(["HIGH", "LOW", "MODIFIER"], list(actual["SNPEFF_TOP_EFFECT_IMPACT"].values))
        self.assertEquals(["1", ".", "0"], list(actual["JQ_SUMMARY_SOM_COUNT|P1|TUMOR"].values))
        self.assertEquals(["1", "1", "."], list(actual["JQ_SUMMARY_SOM_COUNT|P2|TUMOR"].values))

    def test_summarize_dataMatrix(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name = 'SNPEFF_TOP_EFFECT_IMPACT',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        snpEff = rollup.Effect([formatRule], args)

        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\tHIGH\t1\t1
BRCA1\tLOW\t.\t1
CREBBP\tHIGH\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        (data_df, style_dfs) = snpEff.summarize(input_df)

        self.assertEquals(1, len(style_dfs))
        self.assertEquals(data_df.shape, style_dfs[0].shape)

        data_df = data_df.applymap(str)
        expected_header = \
'''GENE_SYMBOL
SnpEff|overall impact rank
SnpEff|overall impact score
SnpEff|impact category|HIGH
SnpEff|impact category|MODERATE
SnpEff|impact category|LOW
SnpEff|impact category|MODIFIER
JQ_SUMMARY_SOM_COUNT|P1|TUMOR
JQ_SUMMARY_SOM_COUNT|P2|TUMOR'''

        expected_string = expected_header.replace("\n", "\t") + '''
BRCA1|1|2000001000|2|0|1|0|h|hl
CREBBP|2|0|0|0|0|0||'''.replace("|", "\t")
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in data_df.values])

    def test_summarize_dataMatrixIgnoresNullOrZeroImpact(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 1)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name = 'SNPEFF_TOP_EFFECT_IMPACT',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        snpEff = rollup.Effect([formatRule], args)

        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\t.\t1\t1
BRCA1\t\t.\t1
CREBBP\tHIGH\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        (data_df, style_dfs) = snpEff.summarize(input_df)

        self.assertEquals(1, len(style_dfs))
        self.assertEquals(data_df.shape, style_dfs[0].shape)

        data_df = data_df.applymap(str)
        expected_header = \
'''GENE_SYMBOL
SnpEff|overall impact rank
SnpEff|overall impact score
SnpEff|impact category|HIGH
SnpEff|impact category|MODERATE
SnpEff|impact category|LOW
SnpEff|impact category|MODIFIER
JQ_SUMMARY_SOM_COUNT|P1|TUMOR
JQ_SUMMARY_SOM_COUNT|P2|TUMOR'''

        expected_string = expected_header.replace("\n", "\t") + '''
CREBBP|1|0|0|0|0|0||'''.replace("|", "\t")
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in data_df.values])

    def test_summarize_formatMatrix(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        formatRule = MockFormatRule(FORMAT_DF)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name = 'SNPEFF_TOP_EFFECT_IMPACT',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        snpEff = rollup.Effect([formatRule], args)

        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\tHIGH\t1\t1
BRCA1\tLOW\t.\t1
CREBBP\tHIGH\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        (data_df, format_dfs) = snpEff.summarize(input_df)

        self.assertIs(data_df, formatRule.last_data_df)
        self.assertEquals(1, len(format_dfs))
        self.assertIs(FORMAT_DF, format_dfs[0])

    def test_summarize_multipleFormatRules(self):
        FORMAT_DF1 = pd.DataFrame([[42] * 8] * 2)
        formatRule1 = MockFormatRule(FORMAT_DF1)

        FORMAT_DF2 = pd.DataFrame([["A"] * 8] * 2)
        formatRule2 = MockFormatRule(FORMAT_DF2)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name = 'SNPEFF_TOP_EFFECT_IMPACT',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        snpEff = rollup.Effect([formatRule1, formatRule2], args)

        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\tHIGH\t1\t1
BRCA1\tLOW\t.\t1
CREBBP\tHIGH\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        (data_df, format_dfs) = snpEff.summarize(input_df)

        self.assertIs(data_df, formatRule1.last_data_df)
        self.assertEquals(2, len(format_dfs))
        self.assertIs(FORMAT_DF1, format_dfs[0])
        self.assertIs(FORMAT_DF2, format_dfs[1])

    def test_calculate_score(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\tLOW\t.\t1
BRCA1\tHIGH\t1\t1
CREBBP\tMODERATE\t0\t.'''
        input_df = dataframe(input_string, sep="\t")

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name='SNPEFF_TOP_EFFECT_IMPACT',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        SnpEff = rollup.Effect([MockFormatRule], args)

        scored_df = SnpEff._calculate_score(input_df)
        self.assertEquals(["BRCA1", "CREBBP"], list(scored_df.index.values))

        self.assertEquals(["h", ""], list(scored_df["JQ_SUMMARY_SOM_COUNT|P1|TUMOR"]))
        self.assertEquals(["hl", ""], list(scored_df["JQ_SUMMARY_SOM_COUNT|P2|TUMOR"]))
        int_scores = [int(i) for i in list(scored_df["effect_annotation|overall_effect_score"].values)]
        self.assertEquals([2000001000, 0], int_scores)

    def test_get_impact_category_counts(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR
BRCA1\thh\tl
CREBBP\tm\t.'''
        input_df = dataframe(input_string, sep="\t")

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name = 'SNPEFF_TOP_EFFECT_IMPACT',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        effect = rollup.Effect([MockFormatRule], args)

        grouped_df = input_df.groupby("GENE_SYMBOL").sum()
        category_df = effect._get_effect_category_counts(grouped_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(category_df.index.values))
        self.assertEquals([2, 0], list(category_df["effect_annotation|effect_category|HIGH"]))
        self.assertEquals([0, 1], list(category_df["effect_annotation|effect_category|MODERATE"]))
        self.assertEquals([1, 0], list(category_df["effect_annotation|effect_category|LOW"]))
        self.assertEquals([0, 0], list(category_df["effect_annotation|effect_category|MODIFIER"]))

    def test_calculate_rank(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT_A|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\teffect_annotation|overall_effect_score
BRCA1\thh\tl\t2000001000
CREBBP\tm\t.\t1000'''
        input_df = dataframe(input_string, sep="\t")

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name='foo',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        SnpEff = rollup.Effect([MockFormatRule], args)

        grouped_df = input_df.groupby("GENE_SYMBOL").sum()
        ranked_df = SnpEff._calculate_rank(grouped_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals(["1", "2"], list(ranked_df["effect_annotation|overall_effect_rank"].values))

    def test_calculate_rank_sortsNumerically(self):
        input_string =\
'''GENE_SYMBOL\teffect_annotation|overall_effect_score
geneA\t1
geneB\t10
geneC\t11
geneD\t2
geneE\t20'''
        input_df = dataframe(input_string, sep="\t")
        input_df = input_df.applymap(str)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name='foo',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        SnpEff = rollup.Effect([MockFormatRule], args)
        ranked_df = SnpEff._calculate_rank(input_df)

        self.assertEquals(["5", "3", "2", "4", "1"], list(ranked_df["effect_annotation|overall_effect_rank"].values))

    def test_calculate_rank_tie(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT_A|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\teffect_annotation|overall_effect_score
CREBBP\tm\t.\t1
TANK\tm\t.\t3
EGFR\tm\t.\t3
BRCA2\tm\t.\t3
BRCA1\t.\tm\t12'''
        input_df = dataframe(input_string, sep="\t")

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name='foo',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        SnpEff = rollup.Effect([MockFormatRule], args)

        grouped_df = input_df.groupby("GENE_SYMBOL").sum()
        ranked_df = SnpEff._calculate_rank(grouped_df)

        self.assertEquals(["BRCA1", "BRCA2", "CREBBP", "EGFR", "TANK"], list(ranked_df.index.values))
        self.assertEquals(["1", "2", "5", "2", "2"], list(ranked_df["effect_annotation|overall_effect_rank"].values))

    def test_change_col_order(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tJQ_SUMMARY_SOM_COUNT|P2|TUMOR\teffect_annotation|overall_effect_score\teffect_annotation|overall_effect_rank
BRCA1\th\thl\t20000\t1
CREBBP\tm\t.\t1\t2'''
        input_df = dataframe(input_string, sep="\t")

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name='foo',
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=False)

        effect = rollup.Effect([MockFormatRule], args)

        indexed_df = input_df.set_index(["GENE_SYMBOL"])
        rearranged_df = effect._change_col_order(indexed_df)
        self.assertEquals(["effect_annotation|overall_effect_rank", "effect_annotation|overall_effect_score", "JQ_SUMMARY_SOM_COUNT|P1|TUMOR", "JQ_SUMMARY_SOM_COUNT|P2|TUMOR"],
                            list(rearranged_df.columns.values))


class SummaryColumnsTestCase(unittest.TestCase):
    def setUp(self):
        rollup._SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
        rollup._GENE_SYMBOL = "GENE_SYMBOL"
        rollup._XLSX = False

    def tearDown(self):
        pass

    def test_summarize(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t1\t1
BRCA1\t.\t1
CREBBP\t0\t.'''
        input_df = dataframe(input_string, sep="\t")

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(.*)\|TUMOR',
                         gene_column_name=rollup._GENE_SYMBOL,
                         input_gene_column_name=rollup._GENE_SYMBOL,
                         tsv=False)

        SummaryColumns = rollup.SummaryColumns(args)
        data_df, style_dfs = SummaryColumns.summarize(input_df)

        expected_string =\
'''GENE_SYMBOL\ttotal_affected_samples\tdistinct_variant_loci\ttotal_mutations
BRCA1\t1\t2\t2.0
CREBBP\t0\t0\t0'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df = expected_df.applymap(str)

        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in data_df.values])
        self.assertEquals(1, len(style_dfs))
        self.assertEquals([["", "", ""], ["", "", ""]], [list(i) for i in style_dfs[0].values])

    def test_calculate_affected_samples(self):
        input_string =\
'''gene_symbol|SampleA|SampleB
BRCA1|1|1
BRCA1|.|1
CREBBP|0|.'''
        input_df = dataframe(input_string)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='(Sample.*)',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        summary_cols = rollup.SummaryColumns(args)
        total_variants = summary_cols.calculate_affected_samples(input_df)

        self.assertEquals([2, 0], list(total_variants))

    def test_calculate_total_mutations_simpleInputMatrix(self):
        input_string =\
'''gene_symbol|SampleA|SampleB
BRCA1|1|1
BRCA1|.|1
CREBBP|0|.'''
        input_df = dataframe(input_string)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='(Sample.*)',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        summary_cols = rollup.SummaryColumns(args)
        total_variants =summary_cols.calculate_total_mutations(input_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(total_variants.index.values))
        self.assertEquals([3, 0], list(total_variants))

    def test_calculate_total_mutations_integerInputMatrix(self):
        input_string =\
'''gene_symbol|SampleA|SampleB
BRCA1|2|1
BRCA1|.|3
CREBBP|0|.'''
        input_df = dataframe(input_string)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='(Sample.*)',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        summary_cols = rollup.SummaryColumns(args)
        total_variants =summary_cols.calculate_total_mutations(input_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(total_variants.index.values))
        self.assertEquals([3, 0], list(total_variants))

    def test_calculate_total_mutations_floatInputMatrix(self):
        input_string =\
'''gene_symbol|SampleA|SampleB
BRCA1|0.5|0.75
BRCA1|.|0.1
CREBBP|0|.'''
        input_df = dataframe(input_string)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='(Sample.*)',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        summary_cols = rollup.SummaryColumns(args)
        total_variants =summary_cols.calculate_total_mutations(input_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(total_variants.index.values))
        self.assertEquals([3, 0], list(total_variants))


    def test_calculate_total_loci(self):
        input_string =\
'''gene_symbol|SampleA|SampleB
BRCA1|1|1
BRCA1|.|1
CREBBP|0|.'''
        input_df = dataframe(input_string)

        args = Namespace(input_file="",
                         output_file="",
                         sample_genotype_column_regex='(Sample.*)',
                         gene_column_name='GENE_SYMBOL',
                         tsv=False)

        summary_cols = rollup.SummaryColumns(args)
        total_loci = summary_cols.calculate_variant_loci(input_df)
        self.assertEquals(["BRCA1", "CREBBP"], list(total_loci.index.values))
        self.assertEquals([2, 0], list(total_loci))


class EffectFormatRuleTestCase(unittest.TestCase):
    def test_format_invalidDataMatrix(self):
        input_string =\
'''GENE_SYMBOL|PATIENT_A|SnpEff_overall_impact_rank
HIGH1|z|1
HIGH2|1|1'''
        data_df = dataframe(input_string)
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup.EffectFormatRule()
        actual_df = rule.format(data_df)

        actual_values = list(actual_df.values)
        self.assertEquals(2, len(actual_values))
        self.assertEquals(2, len(actual_values[0]))
        self.assertEquals(2, len(actual_values[1]))

        pd.isnull(actual_values[0][0])
        pd.isnull(actual_values[0][1])
        pd.isnull(actual_values[1][0])
        pd.isnull(actual_values[1][1])

    def test_format_style(self):
        input_string =\
'''GENE_SYMBOL\tPATIENT_A\tSnpEff_overall_impact_rank
MOD\tmml\t1
HIGH\thhmlx\t1
NULL1\t\t'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])

        rule = rollup.EffectFormatRule()
        actual_df = rule.format(data_df)

        expected_index = ["MOD", "HIGH", "NULL1"]
        self.assertEquals(expected_index, list(actual_df.index))

        expected_patient_cells = [{"font_size": "4", "bg_color": "#3377B9", "font_color": "#3377B9"},
                                  {"font_size": "4", "bg_color": "#003366", "font_color": "#003366"},
                                  numpy.nan]
        self.assertEquals(expected_patient_cells, list(actual_df["PATIENT_A"].values))

        actual_values = list(actual_df["SnpEff_overall_impact_rank"].values)
        self.assertEquals(3, len(actual_values))
        pd.isnull(actual_values[0])
        pd.isnull(actual_values[1])
        pd.isnull(actual_values[2])

    def test_format_standalone(self):
        input_string =\
'''GENE_SYMBOL|PATIENT_A|SnpEff_overall_impact_rank
HIGH1|h|1
HIGH2|hhmlx|1
MOD1|m|1
MOD2|mmlx|1
LOW1|l|1
LOW2|llx|1
MOD1|x|1
MOD2|xx|1
NULL1|.|1'''
        data_df = dataframe(input_string)
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup.EffectFormatRule()
        rule._style = lambda x: x
        actual_df = rule.format(data_df)

        expected_string = \
'''GENE_SYMBOL|PATIENT_A|SnpEff_overall_impact_rank
HIGH1|h|
HIGH2|h|
MOD1|m|
MOD2|m|
LOW1|l|
LOW2|l|
MOD1|x|
MOD2|x|
NULL1||'''
        expected_df = dataframe(expected_string)
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df.fillna("", inplace=True)
        self.assertEquals(list([list(i) for i in expected_df.values]),
                          list([list(i) for i in actual_df.values]))

    def test_style_high(self):
        cell_value = "h"
        rule = rollup.EffectFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "4",
                          "bg_color": "#003366",
                          "font_color": "#003366"}

        self.assertEquals(expected_style, actual_style)

    def test_style_moderate(self):
        cell_value = "m"
        rule = rollup.EffectFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "4",
                          "bg_color": "#3377B9",
                          "font_color": "#3377B9"}

        self.assertEquals(expected_style, actual_style)

    def test_style_low(self):
        cell_value = "l"
        rule = rollup.EffectFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "4",
                          "bg_color": "#91C4E8",
                          "font_color": "#91C4E8"}

        self.assertEquals(expected_style, actual_style)

    def test_style_modifier(self):
        cell_value = "x"
        rule = rollup.EffectFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "4",
                          "bg_color": "#CCCCCC",
                          "font_color": "#CCCCCC"}

        self.assertEquals(expected_style, actual_style)


class dbNSFPFormatRuleTestCase(unittest.TestCase):
    def test_format_style(self):
        input_string =\
'''GENE_SYMBOL\tPATIENT_A\tdbNSFP_annotation|overall_damaging_rank\tdbNSFP_annotation|damaging_total
GENE1\t\t\t1
GENE2\t53\t1\t1
GENE3\t94\t1\t2'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup.dbNSFPFormatRule()
        actual_df = rule.format(data_df)

        expected_index = ["GENE1", "GENE2", "GENE3"]
        self.assertEquals(expected_index, list(actual_df.index))

        expected_patient_cells = [numpy.nan,
                                  {"font_size": "12", "bg_color": "white", "font_color": "#000000"},
                                  {"font_size": "12", "bg_color": "orange", "font_color": "#000000"}]
        self.assertEquals(expected_patient_cells, list(actual_df["PATIENT_A"].values))

        actual_values = list(actual_df["dbNSFP_annotation|overall_damaging_rank"].values)
        self.assertEquals(3, len(actual_values))

        pd.isnull(actual_values[0])
        pd.isnull(actual_values[1])
        pd.isnull(actual_values[2])

    def test_format_standalone(self):
        input_string =\
'''GENE_SYMBOL\tPATIENT_A\tdbNSFP_annotation|overall_damaging_rank\tdbNSFP_annotation|damaging_total
GENE1\t3\t1\t1
GENE2\t53\t1\t1
GENE3\t94\t1\t1
GENE4\t157\t1\t1
GENE5\t33\t1\t1
GENE6\t2\t1\t1'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup.dbNSFPFormatRule()
        rule._style = lambda x: x

        rollup._SAMPLENAME_REGEX = "PATIENT.*"
        actual_df = rule.format(data_df)
        actual_df = actual_df.drop("dbNSFP_annotation|damaging_total", 1)

        expected_string = \
'''GENE_SYMBOL\tPATIENT_A\tdbNSFP_annotation|overall_damaging_rank
GENE1\t\t
GENE2\t\t
GENE3\t\t
GENE4\t\t
GENE5\t\t
GENE6\t\t'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df["PATIENT_A"] = [0, 32, 59, 100, 20, 0]

        self.assertEquals(list([list(i) for i in str(expected_df.values)]),
                          list([list(i) for i in str(actual_df.values)]))

    def test_style_lightest(self):
        cell_value = 0
        rule = rollup.dbNSFPFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "white",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_lighter(self):
        cell_value = 12
        rule = rollup.dbNSFPFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#f2eeee",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_darker(self):
        cell_value = 78
        rule = rollup.dbNSFPFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#e99c4e",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_darkest(self):
        cell_value = 100
        rule = rollup.dbNSFPFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "orange",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)


class RankFormatRuleTestCase(unittest.TestCase):
    def test_format_style(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_annotation|overall_damaging_rank
GENE1\t1
GENE2\t
GENE3\t6'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup.RankFormatRule()
        actual_df = rule.format(data_df)

        expected_index = ["GENE1", "GENE2", "GENE3"]
        self.assertEquals(expected_index, list(actual_df.index))

        expected_patient_cells = [{"font_size": "12", "bg_color": "#e27171", "font_color": "#000000"},
                                  numpy.nan,
                                  {"font_size": "12", "bg_color": "white", "font_color": "#000000"}]
        self.assertEquals(expected_patient_cells, list(actual_df["dbNSFP_annotation|overall_damaging_rank"].values))

    def test_format_styleAllTies(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP|overall damaging rank
GENE1\t1
GENE2\t1
GENE3\t1'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup.RankFormatRule()
        actual_df = rule.format(data_df)

        expected_index = ["GENE1", "GENE2", "GENE3"]
        self.assertEquals(expected_index, list(actual_df.index))

        expected_patient_cells = [numpy.nan,
                                  numpy.nan,
                                  numpy.nan,]
        self.assertEquals(expected_patient_cells, list(actual_df["dbNSFP|overall damaging rank"].values))


    def test_format_standalone(self):
        input_string =\
'''GENE_SYMBOL\teffect_annotation|overall_effect_rank
GENE1\t3
GENE2\t157
GENE3\t52'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup.RankFormatRule()
        rule._style = lambda x: x

        actual_df = rule.format(data_df)

        expected_string = \
'''GENE_SYMBOL\teffect_annotation|overall_effect_rank
GENE1\t0
GENE2\t100
GENE3\t31'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])

        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals(list([list(i) for i in expected_df.values]),
                          list([list(i) for i in actual_df.values]))

    def test_format_rankColumnOnly(self):
        input_string =\
'''GENE_SYMBOL\teffect_annotation|overall_effect_rank\tfoo
GENE1\t3\t1
GENE2\t157\t1
GENE3\t52\t1'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup.RankFormatRule()
        rule._style = lambda x: x

        actual_df = rule.format(data_df)

        expected_string = \
'''GENE_SYMBOL\teffect_annotation|overall_effect_rank\tfoo
GENE1\t0\t
GENE2\t100\t
GENE3\t31\t'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])

        expected_df.fillna("", inplace=True)
        actual_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals(list([list(i) for i in expected_df.values]),
                          list([list(i) for i in actual_df.values]))

    def test_format_multipleRankColumns(self):
        input_string =\
'''GENE_SYMBOL\teffect_annotation|overall_effect_rank\tdbNSFP_annotation|overall_damaging_rank
GENE1\t3\t1
GENE2\t157\t58
GENE3\t52\t93'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup.RankFormatRule()
        rule._style = lambda x: x

        actual_df = rule.format(data_df)

        expected_string = \
'''GENE_SYMBOL\teffect_annotation|overall_effect_rank\tdbNSFP_annotation|overall_damaging_rank
GENE1\t0\t0
GENE2\t100\t61
GENE3\t31\t100'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])

        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals(list([list(i) for i in expected_df.values]),
                          list([list(i) for i in actual_df.values]))

    def test_style_lightest(self):
        cell_value = "100"
        rule = rollup.RankFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "white",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_lighter(self):
        cell_value = "78"
        rule = rollup.RankFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#e9dddd",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_darker(self):
        cell_value = "12"
        rule = rollup.RankFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#f22c2c",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_darkest(self):
        cell_value = "1"
        rule = rollup.RankFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#fe0404",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)


class GeneRollupFunctionalTestCase(unittest.TestCase):
    def setUp(self):
        rollup._SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
        rollup._GENE_SYMBOL = "GENE_SYMBOL"
        rollup._XLSX = False
        rollup._DBNSFP_COLUMN = 'dbNSFP_rollup_damaging'
        rollup._EFFECT_COLUMN = 'SNPEFF_TOP_EFFECT_IMPACT'

    def tearDown(self):
        pass

    def test_rollup(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))

            module_testdir = os.path.join(test_dir, "functional_tests", "gene_rollup")
            input_file = os.path.join(module_testdir, "input", "input.csv")
            output_file = os.path.join(output_dir.path, "rollup.csv")

            expected_file = os.path.join(module_testdir, "benchmark", "rollup.csv")

            dbnsfp_column = 'dbNSFP_rollup_damaging'
            effect_column = 'SNPEFF_TOP_EFFECT_IMPACT'
            rollup._DESIRED_ANNOTATIONS = set([dbnsfp_column, effect_column])

            args = Namespace(input_file=input_file,
                         output_file=output_file,
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(P.)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         dbnsfp_column_name=dbnsfp_column,
                         effect_column_name=effect_column,
                         effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(','),
                         tsv=True)

            rollup._rollup(args)

            expected = open(expected_file).readlines()

            print(expected)
            print(open(output_file).readlines())
            for i, actual in enumerate(open(output_file).readlines()):
                self.maxDiff=None
                self.assertEquals(expected[i], actual)

    def test_rollup_onlySnpEff(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))

            module_testdir = os.path.join(test_dir, "functional_tests", "gene_rollup")
            input_file = os.path.join(module_testdir, "input", "SnpEff_only_input.csv")
            output_file = os.path.join(output_dir.path, "rollup.csv")

            expected_file = os.path.join(module_testdir, "benchmark", "SnpEff_only_rollup.csv")

            effect_column = 'SNPEFF_TOP_EFFECT_IMPACT'
            effect_column_values='HIGH,MODERATE,LOW,MODIFIER'.split(',')
            rollup._DESIRED_ANNOTATIONS = set([effect_column])

            args = Namespace(input_file=input_file,
                         output_file=output_file,
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(P.)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         effect_column_name=effect_column,
                         effect_column_values=effect_column_values,
                         tsv=True)

            rollup._rollup(args)

            expected = open(expected_file).readlines()

            print(expected)
            print(open(output_file).readlines())
            for i, actual in enumerate(open(output_file).readlines()):
                self.assertEquals(expected[i], actual)

    def test_rollup_onlydbNSFP(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))

            module_testdir = os.path.join(test_dir, "functional_tests", "gene_rollup")
            input_file = os.path.join(module_testdir, "input", "dbNSFP_only_input.csv")
            output_file = os.path.join(output_dir.path, "rollup.csv")

            expected_file = os.path.join(module_testdir, "benchmark", "dbNSFP_only_rollup.csv")

            dbnsfp_column = 'dbNSFP_rollup_damaging'
            rollup._DESIRED_ANNOTATIONS = set([dbnsfp_column])

            args = Namespace(input_file=input_file,
                         output_file=output_file,
                         sample_genotype_column_regex='JQ_SUMMARY_SOM_COUNT\|(P.)\|TUMOR',
                         gene_column_name='GENE_SYMBOL',
                         dbnsfp_column_name=dbnsfp_column,
                         tsv=True)

            rollup._rollup(args)

            expected = open(expected_file).readlines()

            print(expected)
            print(open(output_file).readlines())
            for i, actual in enumerate(open(output_file).readlines()):
                self.assertEquals(expected[i], actual)
