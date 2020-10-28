from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from apollo import __version__, _program

setup(name='apollo',
      version=__version__,
      packages=find_packages(),
      scripts=["apollo/scripts/Snakefile",
            "apollo/scripts/paramether.py",
            "apollo/scripts/custom_logger.py",
            "apollo/scripts/log_handler_handle.py",
            "apollo/scripts/apollofunks.py",
            "apollo/scripts/count_cpgs.smk"],
      package_data={"apollo":["data/*",
                  "data/phalacrocorax/*",
                  "data/mus/*",
                  "data/apodemus/*"]},
      install_requires=[
            "biopython>=1.70",
            "numpy>=1.13.3",
            "parasail"
        ],
      description='Predicted Epigenetic Age Clock',
      url='github.com/cov-lineages/apollo',
      author='Aine OToole & Tom Little',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = apollo.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
