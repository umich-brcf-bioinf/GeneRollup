GeneRollup
==========

Command-line tool that accepts a TSV file of variants by samples and emits an
XLSX file of genes by samples, rolling up one or more variants into a single
gene row.

[![Build Status](https://travis-ci.org/umich-brcf-bioinf/GeneRollup.svg?branch=master)](https://travis-ci.org/umich-brcf-bioinf/GeneRollup)

The official repository is at:

https://github.com/umich-brcf-bioinf/GeneRollup

----------
Quickstart
----------
For help and a description of command-line options, run:
   `$ generollup -h`

Rollup variants to genes
========================
There are some constraints on the input file:
- It is a tab-separated file.
- Each row represents a single variant locus associated with a gene. (Recommend
  that you filter the input to exclude intergenic regions.)
- Each row represents a single annotation, typically the "worst" impact.
  Since each locus could be associated with more than one alt allele and/or multiple annotations, you may need to process the variant calling output to limit to the annotation of choice.
- The column headers must include
  - gene symbol
  - sample genotype columns
  - some combination of variant annotation
    - SnpEff: Impact (High,Moderate,Low,Modifier)
    - dbNSFP: count of annotation sources that designated this variant as damaging
    - VarSeq (GoldenHelix)

Here is toy example input file to illustrate the syntax. This example assumes you have annotated with dbNSFP and SnpEff:

------------------------------------------------------------------------------------------------
| #CHROM | POS     |REF| ALT |GENE_SYMBOL|dbNSFP_damaging|SNPEFF_IMPACT|GT.P1.TUMOR|GT.P2.TUMOR|
|--------|---------|---|-----|-----------|---------------|-------------|-----------|-----------|
| chr1   | 14948   | G |  A  | BRCA1     |       1       |    HIGH     |    0/0    |    0/1    |
| chr1   | 137622  | G |  A  | TANK      |       0       |      .      |     .     |    1/1    |
| chr1   | 1147545 | A | G,C | CREBBP    |       7       |    LOW      |    0/2    |    0/0    |
| chr1   | 1147548 | A |  T  | CREBBP    |       7       |  MODERATE   |     .     |    0/1    |
------------------------------------------------------------------------------------------------

~~~
$ generollup \
       sample_genotype_column_regex='GT\.(.*).TUMOR' \
       gene_column_name='GENE_SYMBOL' \
       effect_column_name='SNPEFF_IMPACT' \
       effect_column_values='HIGH,MODERATE,LOW,MODIFIER' \
       dbnsfp_column_name='dbNSFP_damaging' \
       input.tsv output.xlsx
~~~

See documentation on github for more details on how the input variants are
aggregated.

====

Email bfx-jacquard@umich.edu for support and questions.

UM BRCF Bioinformatics Core
