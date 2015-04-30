#pylint: disable=line-too-long, too-many-public-methods, invalid-name
#pylint: disable=maybe-no-member
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
'''#CHROM\tPOS\tREF\tALT\tGENE_SYMBOL\tdbNSFP_rollup_damaging\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT
1\t2\t3\t4\t5\t6\t7\t8'''
        actual_df = rollup_genes._create_df(StringIO(input_string))
        self.assertEquals(["#CHROM", "POS", "REF", "ALT", "GENE_SYMBOL", "dbNSFP_rollup_damaging", "SNPEFF_TOP_EFFECT_IMPACT", "JQ_SUMMARY_SOM_COUNT"],
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
'''#CHROM\tPOS\tREF\tALT\tGENE_SYMBOL\tdbNSFP_rollup_damaging\tfoo
1\t2\t3\t4\t5\t6\t7
1\t2\t3\t4\t5\t6\t7'''
        self.assertRaisesRegexp(BaseException,
                                    "Input file is missing required headers",
                                    rollup_genes._create_df,
                                    StringIO(input_string))

    def test_remove_unnecessary_columns(self):
        input_string =\
'''#CHROM|POS|REF|ALT|GENE_SYMBOL|dbNSFP_rollup_damaging|JQ_SUMMARY_SOM_COUNT_A|JQ_SUMMARY_SOM_COUNT_B|FOO|BAR
1|2|3|4|5|6'''
        input_df = dataframe(input_string)

        actual_df = rollup_genes._remove_unnecessary_columns(input_df)
        self.assertEquals(["#CHROM", "POS", "REF", "ALT", "GENE_SYMBOL", "dbNSFP_rollup_damaging", "JQ_SUMMARY_SOM_COUNT_A", "JQ_SUMMARY_SOM_COUNT_B"],
                          list(actual_df.columns.values))


class dbNSFPTestCase(unittest.TestCase):
    def test_melt_df(self):
        input_string =\
'''#CHROM|POS|REF|ALT|GENE_SYMBOL|dbNSFP_rollup_damaging|SNPEFF_TOP_EFFECT_IMPACT|JQ_SUMMARY_SOM_COUNT_A|JQ_SUMMARY_SOM_COUNT_B
1|2|3|4|BRCA1|4|.|2|3
1|2|3|4|BRCA1|5|.|3|4'''
        input_df = dataframe(input_string)
        dbNSFP = rollup_genes.dbNSFP()
        actual_df = dbNSFP.melt_df(input_df)

        self.assertEquals([".", "2", "BRCA1", "1", "4", "3", "4", "JQ_SUMMARY_SOM_COUNT_A", "2"], list(actual_df.values[0]))
        self.assertEquals([".", "2", "BRCA1", "1", "4", "3", "5", "JQ_SUMMARY_SOM_COUNT_A", "3"], list(actual_df.values[1]))
        self.assertEquals([".", "2", "BRCA1", "1", "4", "3", "4", "JQ_SUMMARY_SOM_COUNT_B", "3"], list(actual_df.values[2]))
        self.assertEquals([".", "2", "BRCA1", "1", "4", "3", "5", "JQ_SUMMARY_SOM_COUNT_B", "4"], list(actual_df.values[3]))

    def test_melt_df_excludeNullGeneSymbols(self):
        input_string =\
'''#CHROM|POS|REF|ALT|GENE_SYMBOL|dbNSFP_rollup_damaging|SNPEFF_TOP_EFFECT_IMPACT|JQ_SUMMARY_SOM_COUNT_A|JQ_SUMMARY_SOM_COUNT_B
1|2|3|4|.|4|.|2|3
1|2|3|4|BRCA1|5|.|3|4'''
        input_df = dataframe(input_string)
        dbNSFP = rollup_genes.dbNSFP()
        actual_df = dbNSFP.melt_df(input_df)
        print actual_df
        self.assertEquals([".", "2", "BRCA1", "1", "4", "3", "5", "JQ_SUMMARY_SOM_COUNT_A", "3"], list(actual_df.values[0]))
        self.assertEquals([".", "2", "BRCA1", "1", "4", "3", "5", "JQ_SUMMARY_SOM_COUNT_B", "4"], list(actual_df.values[1]))

    def calculate_total_samples(self):
        input_string =\
'''#CHROM|POS|REF|ALT|GENE_SYMBOL|dbNSFP_rollup_damaging|Sample|Sample_Data
1|12|A|G|BRCA1|3|JQ_SUMMARY_SOM_COUNT_A|0
1|23|A|G|BRCA1|3|JQ_SUMMARY_SOM_COUNT_B|.
1|42|A|G|CREBBP|3|JQ_SUMMARY_SOM_COUNT_A|0
1|42|A|G|CREBBP|2|JQ_SUMMARY_SOM_COUNT_B|1'''
        input_df = dataframe(input_string)
        dbNSFP = rollup_genes.dbNSFP()

        total_variants = dbNSFP.calculate_total_samples(input_df)

        self.assertEquals([1, 2], list(total_variants))

    def test_calculate_total_variants(self):
        input_string =\
'''#CHROM|POS|REF|ALT|GENE_SYMBOL|dbNSFP_rollup_damaging|Sample|Sample_Data
1|12|A|G|BRCA1|3|JQ_SUMMARY_SOM_COUNT_A|0
1|23|A|G|BRCA1|3|JQ_SUMMARY_SOM_COUNT_B|.
1|42|A|G|CREBBP|3|JQ_SUMMARY_SOM_COUNT_A|0
1|42|A|G|CREBBP|2|JQ_SUMMARY_SOM_COUNT_B|1'''
        input_df = dataframe(input_string)
        dbNSFP = rollup_genes.dbNSFP()

        total_variants = dbNSFP.calculate_total_variants(input_df)

        self.assertEquals([2, 1], list(total_variants))

    def test_pivot_df(self):
        input_string =\
'''GENE_SYMBOL|dbNSFP_rollup_damaging|Sample|Sample_Data
BRCA1|3|JQ_SUMMARY_SOM_COUNT_A|2
BRCA1|4|JQ_SUMMARY_SOM_COUNT_B|3'''
        input_df = dataframe(input_string)
        dbNSFP = rollup_genes.dbNSFP()
        actual_df = dbNSFP.pivot_df(input_df)

        self.assertEquals([("dbNSFP_rollup_damaging", "JQ_SUMMARY_SOM_COUNT_A"), ("dbNSFP_rollup_damaging", "JQ_SUMMARY_SOM_COUNT_B")],
                              list(actual_df.columns.values))
        self.assertEquals(["3", "4"], list(actual_df.values[0]))

    def test_pivot_df_if0(self):
        input_string =\
'''GENE_SYMBOL|dbNSFP_rollup_damaging|Sample|Sample_Data
BRCA1|0|JQ_SUMMARY_SOM_COUNT_A|2
BRCA1|4|JQ_SUMMARY_SOM_COUNT_B|3'''
        input_df = dataframe(input_string)
        dbNSFP = rollup_genes.dbNSFP()
        actual_df = dbNSFP.pivot_df(input_df)

        self.assertEquals([("dbNSFP_rollup_damaging", "JQ_SUMMARY_SOM_COUNT_A"), ("dbNSFP_rollup_damaging", "JQ_SUMMARY_SOM_COUNT_B")],
                              list(actual_df.columns.values))
        self.assertEquals(["0", "4"], list(actual_df.values[0]))

#TODO: determine how to alter dataframe() to account for pivoted df
    def test_calculate_rank(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tSample\tSample_Data
BRCA1\t3\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\t2
BRCA1\t3\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\t3
CREBBP\t7\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\t2'''
        input_df = dataframe(input_string, sep="\t")
        dbNSFP = rollup_genes.dbNSFP()
        pivoted_df = dbNSFP.pivot_df(input_df)

        ranked_df = dbNSFP.calculate_rank(pivoted_df)

        self.assertEquals("BRCA1", ranked_df.ix[1,:].name)
        self.assertEquals("2", ranked_df.ix[1,'dbNSFP|overall damaging rank'].values[0])
        self.assertEquals("CREBBP", ranked_df.ix[0,:].name)
        self.assertEquals("1", ranked_df.ix[0,'dbNSFP|overall damaging rank'].values[0])

    def test_calculate_rank_tie(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tSample\tSample_Data
CREBBP\t7\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\t2
BRCA1\t7\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\t2'''
        input_df = dataframe(input_string, sep="\t")
        dbNSFP = rollup_genes.dbNSFP()
        pivoted_df = dbNSFP.pivot_df(input_df)

        ranked_df = dbNSFP.calculate_rank(pivoted_df)

        self.assertEquals("BRCA1", ranked_df.ix[0,:].name)
        self.assertEquals("1", ranked_df.ix[0,'dbNSFP|overall damaging rank'].values[0])
        self.assertEquals("CREBBP", ranked_df.ix[1,:].name)
        self.assertEquals("1", ranked_df.ix[1,'dbNSFP|overall damaging rank'].values[0])

#TODO: determine how to alter dataframe() to account for pivoted df
    def test_rename_columns(self):
        input_string =\
'''GENE_SYMBOL\tdbNSFP_rollup_damaging\tSample\tSample_Data\ttotal variants
BRCA1\t3\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\t2\t1
BRCA1\t4\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\t3\t1'''
        input_df = dataframe(input_string, sep="\t")
        dbNSFP = rollup_genes.dbNSFP()
        pivoted_df = dbNSFP.pivot_df(input_df)

        rearranged_df = rollup_genes._rearrange_columns(pivoted_df, dbNSFP)
        self.assertEquals("dbNSFP|damaging votes|P1|NORMAL", rearranged_df.columns[0])
        self.assertEquals("dbNSFP|damaging votes|P1|TUMOR", rearranged_df.columns[1])


class SnpEffTestCase(unittest.TestCase):
    def test_melt_df(self):
        input_string =\
'''#CHROM|POS|REF|ALT|GENE_SYMBOL|dbNSFP_rollup_damaging|SNPEFF_TOP_EFFECT_IMPACT|JQ_SUMMARY_SOM_COUNT_A|JQ_SUMMARY_SOM_COUNT_B|total variants
1|2|3|4|BRCA1|.|HIGH|2|3|.
1|2|3|4|BRCA1|.|LOW|3|4|.'''
        input_df = dataframe(input_string)
        SnpEff = rollup_genes.SnpEff()
        actual_df = SnpEff.melt_df(input_df)

        self.assertEquals(["HIGH", "2", "BRCA1", "1", "4", "3", ".", "JQ_SUMMARY_SOM_COUNT_A", "2"], list(actual_df.values[0]))
        self.assertEquals(["LOW", "2", "BRCA1", "1", "4", "3", ".", "JQ_SUMMARY_SOM_COUNT_A", "3"], list(actual_df.values[1]))
        self.assertEquals(["HIGH", "2", "BRCA1", "1", "4", "3", ".", "JQ_SUMMARY_SOM_COUNT_B", "3"], list(actual_df.values[2]))
        self.assertEquals(["LOW", "2", "BRCA1", "1", "4", "3", ".", "JQ_SUMMARY_SOM_COUNT_B", "4"], list(actual_df.values[3]))

    def test_melt_df_excludeNullGeneSymbols(self):
        input_string =\
'''#CHROM|POS|REF|ALT|GENE_SYMBOL|dbNSFP_rollup_damaging|SNPEFF_TOP_EFFECT_IMPACT|JQ_SUMMARY_SOM_COUNT_A|JQ_SUMMARY_SOM_COUNT_B|total variants
1|2|3|4|.|.|LOW|2|3|.
1|2|3|4|BRCA1|.|HIGH|3|4|.'''
        input_df = dataframe(input_string)
        SnpEff = rollup_genes.SnpEff()
        actual_df = SnpEff.melt_df(input_df)

        self.assertEquals(["HIGH", "2", "BRCA1", "1", "4", "3", ".", "JQ_SUMMARY_SOM_COUNT_A", "3"], list(actual_df.values[0]))
        self.assertEquals(["HIGH", "2", "BRCA1", "1", "4", "3", ".", "JQ_SUMMARY_SOM_COUNT_B", "4"], list(actual_df.values[1]))

    def test_pivot_df(self):
        input_string =\
'''GENE_SYMBOL|dbNSFP_rollup_damaging|SNPEFF_TOP_EFFECT_IMPACT|Sample|Sample_Data
BRCA1|.|MODERATE|JQ_SUMMARY_SOM_COUNT_A|2
BRCA1|.|MODIFIER|JQ_SUMMARY_SOM_COUNT_B|3'''
        input_df = dataframe(input_string)
        SnpEff = rollup_genes.SnpEff()
        actual_df = SnpEff.pivot_df(input_df)

        self.assertEquals(["MODERATE", "MODIFIER"],
                              list(actual_df.columns.values))
        self.assertEquals([1, 0], list(actual_df.values[0]))

    def test_pivot_df_ifNull(self):
        input_string =\
'''GENE_SYMBOL|dbNSFP_rollup_damaging|SNPEFF_TOP_EFFECT_IMPACT|Sample|Sample_Data
BRCA1|.|.|JQ_SUMMARY_SOM_COUNT_A|2
BRCA1|.|HIGH|JQ_SUMMARY_SOM_COUNT_B|3'''
        input_df = dataframe(input_string)
        SnpEff = rollup_genes.SnpEff()
        actual_df = SnpEff.pivot_df(input_df)

        self.assertEquals([".","HIGH"],
                              list(actual_df.columns.values))
        self.assertEquals([1, 0], list(actual_df.values[0]))

    def test_calculate_rank(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tSample\tSample_Data
BRCA1\tHIGH\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\t2
BRCA1\tLOW\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\t3
CREBBP\tMODERATE\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\t2'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()
        pivoted_df = SnpEff.pivot_df(input_df)

        ranked_df = SnpEff.calculate_rank(pivoted_df)

        self.assertEquals("BRCA1", ranked_df.ix[0,:].name)
        self.assertEquals("1", ranked_df.ix[0,'SnpEff|overall impact rank'].values[0])
        self.assertEquals("CREBBP", ranked_df.ix[1,:].name)
        self.assertEquals("2", ranked_df.ix[1,'SnpEff|overall impact rank'].values[0])

    def test_calculate_rank_tie(self):
        input_string =\
'''GENE_SYMBOL\tSNPEFF_TOP_EFFECT_IMPACT\tSample\tSample_Data
CREBBP\tLOW\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\t2
BRCA1\tLOW\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\t2'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()
        pivoted_df = SnpEff.pivot_df(input_df)

        ranked_df = SnpEff.calculate_rank(pivoted_df)

        self.assertEquals("BRCA1", ranked_df.ix[0,:].name)
        self.assertEquals("1", ranked_df.ix[0,'SnpEff|overall impact rank'].values[0])
        self.assertEquals("CREBBP", ranked_df.ix[1,:].name)
        self.assertEquals("1", ranked_df.ix[1,'SnpEff|overall impact rank'].values[0])

    def test_rename_columns(self):
        input_string =\
'''#CHROM\tPOS\tREF\tALT\tGENE_SYMBOL\tdbNSFP_rollup_damaging\tSNPEFF_TOP_EFFECT_IMPACT\tJQ_SUMMARY_SOM_COUNT|P1|NORMAL\tJQ_SUMMARY_SOM_COUNT|P1|TUMOR\ttotal impacted samples
1\t2\t3\t4\tBRCA1\t.\tHIGH\t2\t3\t.
1\t2\t3\t4\tBRCA1\t.\tLOW\t3\t4\t.'''
        input_df = dataframe(input_string, sep="\t")
        SnpEff = rollup_genes.SnpEff()
        melted_df = SnpEff.melt_df(input_df)
        pivoted_df = SnpEff.pivot_df(melted_df)
        ranked_df = SnpEff.calculate_rank(pivoted_df)

        rearranged_df = rollup_genes._rearrange_columns(ranked_df, SnpEff)

        self.assertEquals("SnpEff|impact|P1|NORMAL", rearranged_df.columns[1])
        self.assertEquals("SnpEff|impact|P1|TUMOR", rearranged_df.columns[2])


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

