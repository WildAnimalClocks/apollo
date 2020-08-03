from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from peaclock import __version__, _program

setup(name='peaclock',
      version=__version__,
      packages=find_packages(),
      scripts=["peaclock/scripts/Snakefile",
      "peaclock/scripts/parse_paf.py",
      "peaclock/scripts/process_local_trees.smk",
      "peaclock/scripts/make_report.py",
      "peaclock/scripts/report_template.pmd"],
      package_data={"peaclock":["data/reference.fasta",
                              "data/footer.png"]},
      install_requires=[
            "biopython>=1.70",
            "pytools>=2020.1",
            "pweave>=0.30.3",
            "matplotlib>=3.2.1",
            'pandas>=1.0.1',
            'pysam>=0.15.4',
            "scipy>=1.4.1",
            "numpy>=1.13.3"
        ],
      description='Predicted Epigenetic Age Clock',
      url='github.com/cov-lineages/peaclock',
      author='Aine OToole & Tom Little',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = peaclock.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
