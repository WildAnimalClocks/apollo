from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from peaclock import __version__, _program

setup(name='peaclock',
      version=__version__,
      packages=find_packages(),
      scripts=["peaclock/scripts/Snakefile",
            "peaclock/scripts/paramether.py",
            "peaclock/scripts/custom_logger.py",
            "peaclock/scripts/log_handler_handle.py",
            "peaclock/scripts/peaclockfunks.py",
            "peaclock/scripts/count_cpgs.smk"],
      package_data={"peaclock":["data/*",
                  "data/phalacrocorax/*",
                  "data/mus/*",
                  "data/apodemus/*"]},
      install_requires=[
            "biopython>=1.70",
            "numpy>=1.13.3",
            "parasail"
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
