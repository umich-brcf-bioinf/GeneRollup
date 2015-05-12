#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=maybe-no-member, too-few-public-methods
from __future__ import absolute_import

from argparse import Namespace
from StringIO import StringIO
import filecmp
import os
import unittest

from testfixtures import TempDirectory

import geneRollup.rollup_genes as rollup_genes
import pandas as pd


def dataframe(input_data, sep="|"):
    return pd.read_csv(StringIO(input_data), sep=sep, header=0)

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

class GeneRollupArgsTestCase(unittest.TestCase):
    def setUp(self):
        self.default_sample_regex = rollup_genes._SAMPLENAME_REGEX
        self.default_gene = rollup_genes._GENE_SYMBOL
        self.default_xlsx = rollup_genes._XLSX

    def tearDown(self):
        rollup_genes._SAMPLENAME_REGEX = self.default_sample_regex
        rollup_genes._GENE_SYMBOL = self.default_gene
        rollup_genes._XLSX = self.default_xlsx

    def test_change_global_variables_sample(self):
        args = Namespace(input_file="input",
                         output_file="output",
                         sample_column_regex="Sample.*",
                         gene_column_name=None,
                         xlsx=None)
        rollup_genes._change_global_variables(args)
        self.assertEquals("Sample.*", rollup_genes._SAMPLENAME_REGEX)

    def test_change_global_variables_gene(self):
        args = Namespace(input_file="input",
                         output_file="output",
                         sample_column_regex=None,
                         gene_column_name="GeneSymbol",
                         xlsx=None)
        rollup_genes._change_global_variables(args)
        self.assertEquals("GeneSymbol", rollup_genes._GENE_SYMBOL)

    def test_change_global_variables_xlsx(self):
        args = Namespace(input_file="input",
                         output_file="output",
                         sample_column_regex=None,
                         gene_column_name=None,
                         xlsx=True)
        rollup_genes._change_global_variables(args)
        self.assertEquals(True, rollup_genes._XLSX)

class GeneRollupTestCase(unittest.TestCase):
    def test_create_df(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT
1\t2\t3\t4\t5\t6\t7\t8'''
        actual_df = rollup_genes._create_df(StringIO(input_string))
        self.assertEquals(["GENE_SYMBOL", "dbNSFP_rollup_damaging", "SNPEFF_TOP_EFFECT_IMPACT", "JQ_SUMMARY_SOM_COUNT"],
                            list(actual_df.columns.values))

    def test_create_df_invalid(self):
        input_string =\
'''headerA\theaderB
1\t2'''
        self.assertRaisesRegexp(BaseException,
                                "Input file is missing required headers",
                                rollup_genes._create_df,
                                StringIO(input_string))

    def test_create_df_missingSamples(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tfoo
1\t2\t3
1\t2\t3'''
        self.assertRaisesRegexp(BaseException,
                                    "Input file is missing required headers",
                                    rollup_genes._create_df,
                                    StringIO(input_string))

    def test_sort_by_dbnsfp_rank(self):
        input_string =\
'''gene symbol\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tSnpEff|overall impact rank\tdbNSFP|overall damaging rank
BRCA1\th\t7\t2
EGFR\tm\t4\t3
SON\tm\t5\t1
BRCA2\tm\t5\t1
CREBBP\thhh\t6\t1'''
        input_df = dataframe(input_string, sep="\t")
        sorted_df = rollup_genes._sort_by_dbnsfp_rank(input_df)

        self.assertEquals(["BRCA2", "SON", "CREBBP", "BRCA1", "EGFR"], list(sorted_df["gene symbol"].values))
        self.assertEquals([1, 1, 1, 2, 3], list(sorted_df["dbNSFP|overall damaging rank"].values))

    def test_translate_to_excel(self):
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
            style_df["PATIENT_A_dbNSFP"] = [{"font_size": "12", "bg_color": "#ffa500", "font_color": "#6699FF"},
                                            "",
                                            {"font_size": "12", "bg_color": "white", "font_color": "#003366"}]
            style_df.fillna("", inplace=True)

            writer = pd.ExcelWriter(output_file, engine="xlsxwriter")
            rollup_genes._translate_to_excel(data_df, style_df, writer)

            script_dir = os.path.dirname(os.path.realpath(__file__))
            expected_output = os.path.join(script_dir,
                                           "functional_tests",
                                           "translate_to_excel",
                                           "expected_output.xlsx")
            self.assertTrue(filecmp.dircmp(expected_output, output_file))

class dbNSFPTestCase(unittest.TestCase):
    def setUp(self):
        rollup_genes._SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
        rollup_genes._GENE_SYMBOL = "GENE_SYMBOL"
        rollup_genes._XLSX = False

    def tearDown(self):
        pass

    def test_summarize_dataMatrix(self):
        FORMAT_DF = pd.DataFrame([[42] * 3] * 2)
        formatRule = MockFormatRule(FORMAT_DF)
        dbNSFP = rollup_genes.dbNSFP(formatRule)

        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t2\t1\t1
BRCA1\t5\t.\t1
CREBBP\t3\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        input_df = input_df.applymap(str)
        (data_df, style_df) = dbNSFP.summarize(input_df)
        self.assertEquals(data_df.shape, style_df.shape)

        data_df = data_df.applymap(str)
        expected_string =\
'''GENE_SYMBOL\tdbNSFP|overall damaging rank\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t1\t2\t7
CREBBP\t2\t3\t0'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in data_df.values])

    def test_summarize_formatMatrix(self):
        FORMAT_DF = pd.DataFrame([[42] * 3] * 2)
        formatRule = MockFormatRule(FORMAT_DF)
        dbNSFP = rollup_genes.dbNSFP(formatRule)

        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t2\t1\t1
BRCA1\t5\t.\t1
CREBBP\t3\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        (data_df, format_df) = dbNSFP.summarize(input_df)

        self.assertIs(data_df, formatRule.last_data_df)
        self.assertIs(FORMAT_DF, format_df)

    def test_calculate_rank(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t2\t1\t1
BRCA1\t3\t.\t1
CREBBP\t2\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        formatRule = MockFormatRule(FORMAT_DF)
        dbNSFP = rollup_genes.dbNSFP(formatRule)

        ranked_df = dbNSFP._calculate_rank(input_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 2], list(ranked_df["dbNSFP|overall damaging rank"].values))

    def test_calculate_rank_tie(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t5\t.\t1
CREBBP\t5\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        formatRule = MockFormatRule(FORMAT_DF)
        dbNSFP = rollup_genes.dbNSFP(formatRule)

        ranked_df = dbNSFP._calculate_rank(input_df)
        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 1], list(ranked_df["dbNSFP|overall damaging rank"].values))

    def test_change_col_order(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tdbNSFP|overall damaging rank
BRCA1\th\thl\t1
CREBBP\tm\t.\t2'''
        input_df = dataframe(input_string, sep="\t")
        formatRule = MockFormatRule(FORMAT_DF)
        dbNSFP = rollup_genes.dbNSFP(formatRule)

        indexed_df = input_df.set_index(["GENE_SYMBOL"])
        rearranged_df = dbNSFP._change_col_order(indexed_df)
        self.assertEquals(["dbNSFP|overall damaging rank", "JQ_SUMMARY_SOM_COUNT|P1|NORMAL", "JQ_SUMMARY_SOM_COUNT|P1|TUMOR"],
                            list(rearranged_df.columns.values))

class SnpEffTestCase(unittest.TestCase):
    def setUp(self):
        rollup_genes._SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
        rollup_genes._GENE_SYMBOL = "GENE_SYMBOL"
        rollup_genes._XLSX = False

    def tearDown(self):
        pass

    def test_summarize_dataMatrix(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        formatRule = MockFormatRule(FORMAT_DF)
        snpEff = rollup_genes.SnpEff(formatRule)

        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\tHIGH\t1\t1
BRCA1\tLOW\t.\t1
CREBBP\tHIGH\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        (data_df, style_df) = snpEff.summarize(input_df)
        self.assertEquals(data_df.shape, style_df.shape)

        data_df = data_df.applymap(str)
        expected_string =\
'''GENE_SYMBOL\tSnpEff|overall impact rank\tSnpEff|overall impact score\tSnpEff|impact category|HIGH\tSnpEff|impact category|MODERATE\tSnpEff|impact category|LOW\tSnpEff|impact category|MODIFIER\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t1\t200000\t2\t0\t1\t0\th\thl
CREBBP\t2\t100000\t1\t0\t0\t0\th\t'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in data_df.values])

    def test_summarize_formatMatrix(self):
        FORMAT_DF = pd.DataFrame([[42] * 8] * 2)
        formatRule = MockFormatRule(FORMAT_DF)
        snpEff = rollup_genes.SnpEff(formatRule)

        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\tHIGH\t1\t1
BRCA1\tLOW\t.\t1
CREBBP\tHIGH\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        (data_df, format_df) = snpEff.summarize(input_df)

        self.assertIs(data_df, formatRule.last_data_df)
        self.assertIs(FORMAT_DF, format_df)

    def test_calculate_score(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\tLOW\t.\t1
BRCA1\tHIGH\t1\t1
CREBBP\tMODERATE\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff(MockFormatRule)

        scored_df = SnpEff._calculate_score(input_df)
        self.assertEquals(["BRCA1", "CREBBP"], list(scored_df.index.values))
        int_scores = [int(i) for i in list(scored_df["SnpEff|overall impact score"].values)]
        self.assertEquals([200000, 1], int_scores)

        self.assertEquals(["h", "m"], list(scored_df["JQ_SUMMARY_SOM_COUNT|P1|NORMAL"]))
        self.assertEquals(["hl", ""], list(scored_df["JQ_SUMMARY_SOM_COUNT|P1|TUMOR"]))

    def test_get_impact_category_counts(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\thh\tl
CREBBP\tm\t.'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff(MockFormatRule)

        grouped_df = input_df.groupby("GENE_SYMBOL").sum()
        category_df = SnpEff._get_impact_category_counts(grouped_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(category_df.index.values))
        self.assertEquals([2, 0], list(category_df["SnpEff|impact category|HIGH"]))
        self.assertEquals([0, 1], list(category_df["SnpEff|impact category|MODERATE"]))
        self.assertEquals([1, 0], list(category_df["SnpEff|impact category|LOW"]))
        self.assertEquals([0, 0], list(category_df["SnpEff|impact category|MODIFIER"]))

    def test_calculate_rank(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT_A|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tSnpEff|overall impact score
BRCA1\thh\tl\t200000
CREBBP\tm\t.\t1'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff(MockFormatRule)

        grouped_df = input_df.groupby("GENE_SYMBOL").sum()
        ranked_df = SnpEff._calculate_rank(grouped_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals(["1", "2"], list(ranked_df["SnpEff|overall impact rank"].values))

    def test_calculate_rank_tie(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT_A|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tSnpEff|overall impact score
BRCA1\t.\tm\t1
CREBBP\tm\t.\t1'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff(MockFormatRule)

        grouped_df = input_df.groupby("GENE_SYMBOL").sum()
        ranked_df = SnpEff._calculate_rank(grouped_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals(["1", "1"], list(ranked_df["SnpEff|overall impact rank"].values))

    def test_change_col_order(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tSnpEff|overall impact score\tSnpEff|overall impact rank
BRCA1\th\thl\t20000\t1
CREBBP\tm\t.\t1\t2'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff(MockFormatRule)

        indexed_df = input_df.set_index(["GENE_SYMBOL"])
        rearranged_df = SnpEff._change_col_order(indexed_df)
        self.assertEquals(["SnpEff|overall impact rank", "SnpEff|overall impact score", "JQ_SUMMARY_SOM_COUNT|P1|NORMAL", "JQ_SUMMARY_SOM_COUNT|P1|TUMOR"],
                            list(rearranged_df.columns.values))


class SummaryColumnsTestCase(unittest.TestCase):
    def setUp(self):
        rollup_genes._SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
        rollup_genes._GENE_SYMBOL = "GENE_SYMBOL"
        rollup_genes._XLSX = False

    def tearDown(self):
        pass

    def test_summarize(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t1\t1
BRCA1\t.\t1
CREBBP\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        SummaryColumns = rollup_genes.SummaryColumns()
        data_df, style_df = SummaryColumns.summarize(input_df)

        expected_string =\
'''GENE_SYMBOL\ttotal impacted samples\tdistinct loci\ttotal mutations
BRCA1\t2\t2\t3
CREBBP\t1\t1\t1'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df = expected_df.applymap(str)

        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in data_df.values])
        self.assertEquals([["", "", ""], ["", "", ""]], [list(i) for i in style_df.values])

    def test_calculate_total_samples(self):
        input_string =\
'''GENE_SYMBOL|SampleA|SampleB
BRCA1|1|1
BRCA1|.|1
CREBBP|0|.'''
        input_df = dataframe(input_string)
        summary_cols = rollup_genes.SummaryColumns()
        total_variants = summary_cols.calculate_total_samples(input_df)

        self.assertEquals([2, 1], list(total_variants))

    def test_calculate_total_mutations(self):
        input_string =\
'''GENE_SYMBOL|SampleA|SampleB
BRCA1|1|1
BRCA1|.|1
CREBBP|0|.'''
        input_df = dataframe(input_string)
        summary_cols = rollup_genes.SummaryColumns()
        total_variants =summary_cols.calculate_total_mutations(input_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(total_variants.index.values))
        self.assertEquals([3, 1], list(total_variants))

    def test_calculate_total_loci(self):
        input_string =\
'''GENE_SYMBOL|SampleA|SampleB
BRCA1|1|1
BRCA1|.|1
CREBBP|0|.'''
        input_df = dataframe(input_string)
        summary_cols = rollup_genes.SummaryColumns()
        total_loci = summary_cols.calculate_total_loci(input_df)
        self.assertEquals(["BRCA1", "CREBBP"], list(total_loci.index.values))
        self.assertEquals([2, 1], list(total_loci))


class SnpEffFormatRuleTestCase(unittest.TestCase):
    def test_format_invalidDataMatrix(self):
        input_string =\
'''GENE_SYMBOL|PATIENT_A|SnpEff_overall_impact_rank
HIGH1|z|1
HIGH2|1|1'''
        data_df = dataframe(input_string)
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup_genes.SnpEffFormatRule()
        actual_df = rule.format(data_df)

        expected_string = \
'''GENE_SYMBOL|PATIENT_A|SnpEff_overall_impact_rank
HIGH1||
HIGH2||'''
        expected_df = dataframe(expected_string)
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df.fillna("", inplace=True)
        self.assertEquals(list([list(i) for i in expected_df.values]),
                          list([list(i) for i in actual_df.values]))

    def test_format_style(self):
        input_string =\
'''GENE_SYMBOL|PATIENT_A|SnpEff_overall_impact_rank
MOD|mml|1
HIGH|hhmlx|1
NULL1||'''
        data_df = dataframe(input_string)
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup_genes.SnpEffFormatRule()
        actual_df = rule.format(data_df)

        expected_index = ["MOD", "HIGH", "NULL1"]
        self.assertEquals(expected_index, list(actual_df.index))

        expected_patient_cells = [{"font_size": "4", "bg_color": "#6699FF", "font_color": "#6699FF"},
                                  {"font_size": "4", "bg_color": "#003366", "font_color": "#003366"},
                                  ""]
        self.assertEquals(expected_patient_cells, list(actual_df["PATIENT_A"].values))

        expected_rank_cells = ["", "", ""]
        self.assertEquals(expected_rank_cells, list(actual_df["SnpEff_overall_impact_rank"].values))

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
        rule = rollup_genes.SnpEffFormatRule()
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
        rule = rollup_genes.SnpEffFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "4",
                          "bg_color": "#003366",
                          "font_color": "#003366"}

        self.assertEquals(expected_style, actual_style)

    def test_style_moderate(self):
        cell_value = "m"
        rule = rollup_genes.SnpEffFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "4",
                          "bg_color": "#6699FF",
                          "font_color": "#6699FF"}

        self.assertEquals(expected_style, actual_style)

    def test_style_low(self):
        cell_value = "l"
        rule = rollup_genes.SnpEffFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "4",
                          "bg_color": "#CCCCFF",
                          "font_color": "#CCCCFF"}

        self.assertEquals(expected_style, actual_style)

    def test_style_modifier(self):
        cell_value = "x"
        rule = rollup_genes.SnpEffFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "4",
                          "bg_color": "#999999",
                          "font_color": "#999999"}

        self.assertEquals(expected_style, actual_style)


class dbNSFPFormatRuleTestCase(unittest.TestCase):
    def test_format_style(self):
        input_string =\
'''GENE_SYMBOL\tPATIENT_A\tdbNSFP|overall damaging rank
GENE1\t\t
GENE2\t53\t1
GENE3\t94\t1'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup_genes.dbNSFPFormatRule()
        actual_df = rule.format(data_df)

        expected_index = ["GENE1", "GENE2", "GENE3"]
        self.assertEquals(expected_index, list(actual_df.index))

        expected_patient_cells = ["",
                                  {"font_size": "12", "bg_color": "#fff", "font_color": "#000000"},
                                  {"font_size": "12", "bg_color": "#ffa500", "font_color": "#000000"}]
        self.assertEquals(expected_patient_cells, list(actual_df["PATIENT_A"].values))

        expected_rank_cells = ["", "", ""]
        self.assertEquals(expected_rank_cells, list(actual_df["dbNSFP|overall damaging rank"].values))

    def test_format_standalone(self):
        input_string =\
'''GENE_SYMBOL\tPATIENT_A\tdbNSFP|overall damaging rank
GENE1\t3\t1
GENE2\t53\t1
GENE3\t94\t1
GENE4\t157\t1
GENE5\t33\t1
GENE6\t2\t1'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup_genes.dbNSFPFormatRule()
        rule._style = lambda x: x

        rollup_genes._SAMPLENAME_REGEX = "PATIENT.*"
        actual_df = rule.format(data_df)

        expected_string = \
'''GENE_SYMBOL\tPATIENT_A\tdbNSFP|overall damaging rank
GENE1\t\t
GENE2\t\t
GENE3\t\t
GENE4\t\t
GENE5\t\t
GENE6\t\t'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df["PATIENT_A"] = ["0", "32", "59", "100", "20", "0"]

        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals(list([list(i) for i in expected_df.values]),
                          list([list(i) for i in actual_df.values]))

    def test_style_lightest(self):
        cell_value = 0
        rule = rollup_genes.dbNSFPFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#fff",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_lighter(self):
        cell_value = "12"
        rule = rollup_genes.dbNSFPFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#f2eeee",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_darker(self):
        cell_value = "78"
        rule = rollup_genes.dbNSFPFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#e99c4e",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_darkest(self):
        cell_value = "100"
        rule = rollup_genes.dbNSFPFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#ffa500",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)


class RankFormatRuleTestCase(unittest.TestCase):
    def test_format_style(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP|overall damaging rank
GENE1\t1
GENE2\t
GENE3\t6'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup_genes.RankFormatRule()
        actual_df = rule.format(data_df)

        expected_index = ["GENE1", "GENE2", "GENE3"]
        self.assertEquals(expected_index, list(actual_df.index))

        expected_patient_cells = [{"font_size": "12", "bg_color": "red", "font_color": "#000000"},
                                  "",
                                  {"font_size": "12", "bg_color": "white", "font_color": "#000000"}]
        self.assertEquals(expected_patient_cells, list(actual_df["dbNSFP|overall damaging rank"].values))

    def test_format_standalone(self):
        input_string =\
'''GENE_SYMBOL\tSnpEff|overall impact rank
GENE1\t3
GENE2\t157
GENE3\t52'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup_genes.RankFormatRule()
        rule._style = lambda x: x

        actual_df = rule.format(data_df)

        expected_string = \
'''GENE_SYMBOL\tSnpEff|overall impact rank
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
'''GENE_SYMBOL\tSnpEff|overall impact rank\tfoo
GENE1\t3\t1
GENE2\t157\t1
GENE3\t52\t1'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup_genes.RankFormatRule()
        rule._style = lambda x: x

        actual_df = rule.format(data_df)

        expected_string = \
'''GENE_SYMBOL\tSnpEff|overall impact rank\tfoo
GENE1\t0\t
GENE2\t100\t
GENE3\t31\t'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])

        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals(list([list(i) for i in expected_df.values]),
                          list([list(i) for i in actual_df.values]))

    def test_format_multipleRankColumns(self):
        input_string =\
'''GENE_SYMBOL\tSnpEff|overall impact rank\tdbNSFP|overall damaging rank
GENE1\t3\t1
GENE2\t157\t58
GENE3\t52\t93'''
        data_df = dataframe(input_string, sep="\t")
        data_df = data_df.set_index(["GENE_SYMBOL"])
        rule = rollup_genes.RankFormatRule()
        rule._style = lambda x: x

        actual_df = rule.format(data_df)

        expected_string = \
'''GENE_SYMBOL\tSnpEff|overall impact rank\tdbNSFP|overall damaging rank
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
        rule = rollup_genes.RankFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "white",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_lighter(self):
        cell_value = "78"
        rule = rollup_genes.RankFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#e9dddd",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_darker(self):
        cell_value = "12"
        rule = rollup_genes.RankFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#f22c2c",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)

    def test_style_darkest(self):
        cell_value = "1"
        rule = rollup_genes.RankFormatRule()
        actual_style = rule._style(cell_value)
        expected_style = {"font_size": "12",
                          "bg_color": "#fe0404",
                          "font_color": "#000000"}

        self.assertEquals(expected_style, actual_style)


class GeneRollupFunctionalTestCase(unittest.TestCase):
    def setUp(self):
        rollup_genes._SAMPLENAME_REGEX = "JQ_SUMMARY_SOM_COUNT.*"
        rollup_genes._GENE_SYMBOL = "GENE_SYMBOL"
        rollup_genes._XLSX = False

    def tearDown(self):
        pass

    def test_rollup_genes(self):
        with TempDirectory() as output_dir:
            test_dir = os.path.dirname(os.path.realpath(__file__))

            module_testdir = os.path.join(test_dir, "functional_tests", "gene_rollup")
            input_file = os.path.join(module_testdir, "input", "input.csv")
            output_file = os.path.join(output_dir.path, "rollup.csv")

            expected_file = os.path.join(module_testdir, "benchmark", "rollup.csv")

            rollup_genes._rollup(input_file, output_file)

            expected = open(expected_file).readlines()

            for i, actual in enumerate(open(output_file).readlines()):
                self.assertEquals(expected[i], actual)

