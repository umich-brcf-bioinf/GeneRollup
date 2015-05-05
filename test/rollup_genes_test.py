#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=maybe-no-member, too-few-public-methods
from __future__ import absolute_import

from argparse import Namespace
from StringIO import StringIO
import os
import unittest

from testfixtures import TempDirectory

import geneRollup.rollup_genes as rollup_genes
import pandas as pd


def dataframe(input_data, sep="|"):
    return pd.read_csv(StringIO(input_data), sep=sep, header=0)

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
        self.assertEquals("JQ_SUMMARY_SOM_COUNT.*", rollup_genes._SAMPLENAME_REGEX)
        args = Namespace(input_file="input",
                         output_file="output",
                         sample_column_regex="Sample.*",
                         gene_column_name=None,
                         xlsx=None)
        rollup_genes._change_global_variables(args)
        self.assertEquals("Sample.*", rollup_genes._SAMPLENAME_REGEX)

    def test_change_global_variables_gene(self):
        self.assertEquals("gene symbol", rollup_genes._GENE_SYMBOL)
        args = Namespace(input_file="input",
                         output_file="output",
                         sample_column_regex=None,
                         gene_column_name="GeneSymbol",
                         xlsx=None)
        rollup_genes._change_global_variables(args)
        self.assertEquals("GeneSymbol", rollup_genes._GENE_SYMBOL)

    def test_change_global_variables_xlsx(self):
        self.assertEquals(False, rollup_genes._XLSX)
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
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tdbNSFP|overall damaging rank
BRCA1\th\thl\t2
EGFR\tm\t.\t3
CREBBP\thhh\t.\t1'''
        input_df = dataframe(input_string, sep="\t")
        input_df = input_df.set_index("GENE_SYMBOL")
        sorted_df = rollup_genes._sort_by_dbnsfp_rank(input_df)

        self.assertEquals([1, 2, 3], list(sorted_df["dbNSFP|overall damaging rank"].values))
        self.assertEquals(["CREBBP", "BRCA1", "EGFR"], list(sorted_df.index.values))

class dbNSFPTestCase(unittest.TestCase):
    def test_summarize(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t2\t1\t1
BRCA1\t3\t.\t1
CREBBP\t2\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        dbNSFP = rollup_genes.dbNSFP()
        summarized_df = dbNSFP.summarize(input_df)

        expected_string =\
'''GENE_SYMBOL\tdbNSFP|overall damaging rank\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t1\t5\t5
CREBBP\t2\t2\t2'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])

        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in summarized_df.values])

    def test_calculate_rank(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t2\t1\t1
BRCA1\t3\t.\t1
CREBBP\t2\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        dbNSFP = rollup_genes.dbNSFP()

        ranked_df = dbNSFP._calculate_rank(input_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 2], list(ranked_df["dbNSFP|overall damaging rank"].values))

    def test_calculate_rank_tie(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t2\t1\t1
BRCA1\t3\t.\t1
CREBBP\t5\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        dbNSFP = rollup_genes.dbNSFP()

        ranked_df = dbNSFP._calculate_rank(input_df)
        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 1], list(ranked_df["dbNSFP|overall damaging rank"].values))

    def test_change_col_order(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tdbNSFP|overall damaging rank
BRCA1\th\thl\t1
CREBBP\tm\t.\t2'''
        input_df = dataframe(input_string, sep="\t")
        dbNSFP = rollup_genes.dbNSFP()

        indexed_df = input_df.set_index(["GENE_SYMBOL"])
        rearranged_df = dbNSFP._change_col_order(indexed_df)
        self.assertEquals(["dbNSFP|overall damaging rank", "JQ_SUMMARY_SOM_COUNT|P1|NORMAL", "JQ_SUMMARY_SOM_COUNT|P1|TUMOR"],
                            list(rearranged_df.columns.values))

class SnpEffTestCase(unittest.TestCase):
    def test_summarize(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\tHIGH\t1\t1
BRCA1\tLOW\t.\t1
CREBBP\tHIGH\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()
        summarized_df = SnpEff.summarize(input_df)
        summarized_df = summarized_df.applymap(str)

        expected_string =\
'''GENE_SYMBOL\tSnpEff|overall impact rank\tSnpEff|overall impact score\tSnpEff|impact category|HIGH\tSnpEff|impact category|MODERATE\tSnpEff|impact category|LOW\tSnpEff|impact category|MODIFIER\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t1.0\t200000.00001\t2\t0\t1\t0\th\thl
CREBBP\t2.0\t100000.0\t1\t0\t0\t0\th\t'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])
        expected_df.fillna("", inplace=True)
        expected_df = expected_df.applymap(str)

        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in summarized_df.values])

    def test_calculate_score(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT_A|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\tHIGH\t1\t1
BRCA1\tLOW\t.\t1
CREBBP\tMODERATE\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()

        scored_df = SnpEff._calculate_score(input_df)
        self.assertEquals(["BRCA1", "CREBBP"], list(scored_df.index.values))
        int_scores = [int(i) for i in list(scored_df["SnpEff|overall impact score"].values)]
        self.assertEquals([200000, 1], int_scores)

    def test_get_impact_category_counts(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT_A|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\thh\tl
CREBBP\tm\t.'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()

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
        SnpEff = rollup_genes.SnpEff()

        grouped_df = input_df.groupby("GENE_SYMBOL").sum()
        ranked_df = SnpEff._calculate_rank(grouped_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 2], list(ranked_df["SnpEff|overall impact rank"].values))

    def test_calculate_rank_tie(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT_A|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tSnpEff|overall impact score
BRCA1\t.\tm\t1
CREBBP\tm\t.\t1'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()

        grouped_df = input_df.groupby("GENE_SYMBOL").sum()
        ranked_df = SnpEff._calculate_rank(grouped_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 1], list(ranked_df["SnpEff|overall impact rank"].values))

    def test_change_col_order(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\tSnpEff|overall impact score\tSnpEff|overall impact rank
BRCA1\th\thl\t20000\t1
CREBBP\tm\t.\t1\t2'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()

        indexed_df = input_df.set_index(["GENE_SYMBOL"])
        rearranged_df = SnpEff._change_col_order(indexed_df)
        self.assertEquals(["SnpEff|overall impact rank", "SnpEff|overall impact score", "JQ_SUMMARY_SOM_COUNT|P1|NORMAL", "JQ_SUMMARY_SOM_COUNT|P1|TUMOR"],
                            list(rearranged_df.columns.values))

class SummaryColumnsTestCase(unittest.TestCase):
    def test_summarize(self):
        input_string =\
'''GENE_SYMBOL\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\t1\t1
BRCA1\t.\t1
CREBBP\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        SummaryColumns = rollup_genes.SummaryColumns()
        summarized_df = SummaryColumns.summarize(input_df)

        expected_string =\
'''GENE_SYMBOL\ttotal impacted samples\tdistinct loci\ttotal mutations
BRCA1\t2\t2\t3
CREBBP\t1\t1\t1'''
        expected_df = dataframe(expected_string, sep="\t")
        expected_df = expected_df.set_index(["GENE_SYMBOL"])

        self.assertEquals([list(i) for i in expected_df.values], [list(i) for i in summarized_df.values])

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


class GeneRollupFunctionalTestCase(unittest.TestCase):
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

