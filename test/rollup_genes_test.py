from __future__ import absolute_import

from argparse import Namespace
import geneRollup.rollup_genes as rollup_genes
import os
from testfixtures import TempDirectory
import unittest


class GeneRollupTestCase(unittest.TestCase):
    def test_create_df(self):
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "headerA\theaderB\n1\t2")
            input_file = os.path.join(input_dir.path, "input.tsv")
            actual_df = rollup_genes._create_df(input_file)

            self.assertEquals(["headerA", "headerB"], list(actual_df.columns.values))

    def test_validate_df_invalid(self):
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "headerA\theaderB\n1\t2")
            input_file = os.path.join(input_dir.path, "input.tsv")
            actual_df = rollup_genes._create_df(input_file)

            self.assertRaisesRegexp(BaseException,
                                    "Input file is missing required headers",
                                    rollup_genes._validate_df,
                                    actual_df)

    def test_validate_df_missingSamples(self):
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "GENE_SYMBOL\tdbNSFP_rollup_damaging\tfoo\n1\t2")
            input_file = os.path.join(input_dir.path, "input.tsv")
            actual_df = rollup_genes._create_df(input_file)

            self.assertRaisesRegexp(BaseException,
                                    "Input file is missing required headers",
                                    rollup_genes._validate_df,
                                    actual_df)

    def test_validate_df_valid(self):
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_CONS_SOM\n1\t2")
            input_file = os.path.join(input_dir.path, "input.tsv")
            actual_df = rollup_genes._create_df(input_file)

            rollup_genes._validate_df(actual_df)
            self.assertTrue(1)

    def test_remove_unnecessary_columns(self):
        with TempDirectory() as input_dir:
            input_dir.write("input.tsv", "GENE_SYMBOL\tdbNSFP_rollup_damaging\tJQ_CONS_SOM_A\tJQ_CONS_SOM_B\tFOO\tBAR\n1\t2")
            input_file = os.path.join(input_dir.path, "input.tsv")
            df = rollup_genes._create_df(input_file)

            actual_df = rollup_genes._remove_unnecessary_columns(df)
            self.assertEquals(["GENE_SYMBOL", "dbNSFP_rollup_damaging", "JQ_CONS_SOM_A", "JQ_CONS_SOM_B"],
                              list(actual_df.columns.values))

