Installing GeneRollup
=====================
GeneRollup has been tested with Python 3 on \*nix, Windows, and OSX.

Prerequisites
-------------
* pandas
* colour
* nosetests, testfixtures (3.0.2) are required for running
  automated tests
* Note that pip installs all required libraries; see [Installing] below.

Installing
----------
You can install from source from github:

``$ pip install git+https://github.com/umich-brcf-bioinf/GeneRollup``

If you don't have root permissions, you can install locally:

``$ pip install git+https://github.com/umich-brcf-bioinf/GeneRollup --user``

Following the pip install, you may need to adjust your path settings to include
home/.local/bin.


If you have already installed the above prerequisites, you can also clone from
github and run directly from the source like so:

``$ git clone https://github.com/umich-brcf-bioinf/GeneRollup``

``$ GeneRollup/rollup-runner.py input.tsv output.xlsx``
