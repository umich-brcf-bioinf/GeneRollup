#pylint: disable=line-too-long, too-many-public-methods, invalid-name
from __future__ import absolute_import

from StringIO import StringIO
import os
import unittest

from testfixtures import TempDirectory

import geneRollup.rollup_genes as rollup_genes
import pandas as pd

def dataframe(input_data, sep="|", index_col=None):
    return pd.read_csv(StringIO(input_data), sep=sep, header=False, dtype='str', index_col=index_col)

class GeneRollupTestCase(unittest.TestCase):
    def test_create_df(self):
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_CONS_SOM\n1\t2")
            input_file = os.path.join(input_dir.path, "input.tsv")
            actual_df = rollup_genes._create_df(input_file)
            #pylint: disable=maybe-no-member
            self.assertEquals(["GENE_SYMBOL", "dbNSFP_rollup_damaging", "JQ_CONS_SOM"],
                              list(actual_df.columns.values))

    def test_create_df_invalid(self):
#         input_string =\
# '''
# headerA|headerB
# 1|2
# '''
#         input_df = dataframe(input_string)
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "headerA\theaderB\n1\t2")
            input_file = os.path.join(input_dir.path, "input.tsv")
            self.assertRaisesRegexp(BaseException,
                                    "Input file is missing required headers",
                                    rollup_genes._create_df,
                                    input_file)

    def test_create_df_missingSamples(self):
#         input_string =\
# '''
# GENE_SYMBOL|dbNSFP_rollup_damaging|foo
# 1|2|3
# '''
#         input_df = dataframe(input_string)
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "GENE_SYMBOL\tdbNSFP_rollup_damaging\tfoo\n1\t2\t3")
            input_file = os.path.join(input_dir.path, "input.tsv")
            self.assertRaisesRegexp(BaseException,
                                    "Input file is missing required headers",
                                    rollup_genes._create_df,
                                    input_file)

    def test_remove_unnecessary_columns(self):
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_CONS_SOM_A\tJQ_CONS_SOM_B\tFOO\tBAR\n1\t2")
            input_file = os.path.join(input_dir.path, "input.tsv")
            df = rollup_genes._create_df(input_file)

            actual_df = rollup_genes._remove_unnecessary_columns(df)
            self.assertEquals(["GENE_SYMBOL", "dbNSFP_rollup_damaging", "JQ_CONS_SOM_A", "JQ_CONS_SOM_B"],
                              list(actual_df.columns.values))

    def test_melt_df(self):
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_CONS_SOM_A\tJQ_CONS_SOM_B\nBRCA1\t4\t2\t3")
            input_file = os.path.join(input_dir.path, "input.tsv")
            df = rollup_genes._create_df(input_file)

            actual_df = rollup_genes._melt_df(df)
            self.assertEquals(["BRCA1", "4", "JQ_CONS_SOM_A", "2"], list(actual_df.values[0]))
            self.assertEquals(["BRCA1", "4", "JQ_CONS_SOM_B", "3"], list(actual_df.values[1]))

    def test_pivot_df(self):
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_CONS_SOM_Sample\tJQ_CONS_SOM_Sample_Data\nBRCA1\t3\tJQ_CONS_SOM_A\t2\nBRCA1\t4\tJQ_CONS_SOM_B\t3")
            input_file = os.path.join(input_dir.path, "input.tsv")
            df = rollup_genes._create_df(input_file)

            actual_df = rollup_genes._pivot_df(df)
            self.assertEquals(["dbNSFP_rollup_damaging"],
                              list(actual_df.columns.values))
            self.assertEquals([3], list(actual_df.values[0]))


