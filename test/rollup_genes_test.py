#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=maybe-no-member, too-few-public-methods
from __future__ import absolute_import

from StringIO import StringIO
import os
import unittest

from testfixtures import TempDirectory

import geneRollup.rollup_genes as rollup_genes
import pandas as pd


def dataframe(input_data, sep="|"):
    return pd.read_csv(StringIO(input_data), sep=sep, header=0)

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

class dbNSFPTestCase(unittest.TestCase):
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
        print ranked_df
        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 1], list(ranked_df["dbNSFP|overall damaging rank"].values))

class SnpEffTestCase(unittest.TestCase):
    def test_calculate_rank(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT_A|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\tHIGH\t1\t1
BRCA1\tLOW\t.\t1
CREBBP\tMODERATE\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()

        ranked_df = SnpEff._calculate_rank(input_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 2], list(ranked_df["SnpEff|overall impact rank"].values))

    def test_calculate_rank_tie(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR
BRCA1\tLOW\t1\t1
CREBBP\tLOW\t0\t.'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()

        ranked_df = SnpEff._calculate_rank(input_df)

        self.assertEquals(["BRCA1", "CREBBP"], list(ranked_df.index.values))
        self.assertEquals([1, 1], list(ranked_df["SnpEff|overall impact rank"].values))

class SummaryColumnsTestCase(unittest.TestCase):
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

