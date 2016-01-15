==========
GeneRollup
==========

Command-line tool that accepts a TSV file of variants by samples and emits an
XLSX file of genes by samples, rolling up one or more variants into a single
gene row.

.. image:: https://travis-ci.org/umich-brcf-bioinf/GeneRollup.svg?branch=develop
    :target: https://travis-ci.org/umich-brcf-bioinf/GeneRollup
    :alt: Build Status

.. image:: https://coveralls.io/repos/umich-brcf-bioinf/GeneRollup/badge.svg?branch=develop&service=github
    :target: https://coveralls.io/github/umich-brcf-bioinf/GeneRollup?branch=develop
    :alt: Coverage Status

The official repository is at:

https://github.com/umich-brcf-bioinf/GeneRollup

----------
Quickstart
----------

Rollup a tab-separated file of variants
=======================================

   $ rollup input.tsv output.xlsx

The input file has to have been annotated by either dbNSFP, SnpEff, or both. By
default, the samples are determined based on a prefix added by Jacquard.


For help and a description of command-line options, run:

   $ rollup -h

====

Email bfx-jacquard@umich.edu for support and questions.

UM BRCF Bioinformatics Core 
