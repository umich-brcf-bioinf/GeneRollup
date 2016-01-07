from setuptools import setup
from setuptools import find_packages
import os
import generollup

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as filename:
        return filename.read()

setup(name='generollup',
      version=generollup.__version__,
      description=('Command-line tool that accepts a TSV file of variants by '
                   'samples and emits an XLSX file of genes by samples, '
                   'rolling up one or more variants into a single gene row.'),
      long_description=(read('README.rst') + '\n\n' +
                        read('CHANGELOG.rst') + '\n\n' +
                        read('AUTHORS.rst')),
      url='https://github.com/umich-brcf-bioinf/GeneRollup',
      author='University of Michigan Bioinformatics Core',
      author_email='bfx-jacquard@umich.edu',
      license='Apache',
      packages=find_packages(exclude=['test*']),
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: Apache Software License',
                   'Operating System :: Unix',
                   'Programming Language :: Python :: 2.7',
                   'Topic :: Scientific/Engineering :: Bio-Informatics'],
      keywords='VCF bioinformatic exome-seq DNA-seq variant-call-format',
      install_requires=['colour', 'pandas'],
      entry_points={'console_scripts': ['generollup=generollup.rollup:main']},
      test_suite='nose.collector',
      tests_require=['nose', 'testfixtures'],
      zip_safe=False)
