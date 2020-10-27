#!/usr/bin/env python3
from peaclock import __version__
import setuptools
import argparse
import os.path
import snakemake
import sys
from tempfile import gettempdir
import tempfile
import pprint
import json
import csv
import os
from datetime import datetime
from Bio import SeqIO
import csv

import pkg_resources
from . import _program

import peaclockfunks as qcfunk
import custom_logger as custom_logger
import log_handler_handle as lh

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description=qcfunk.preamble(__version__), 
    usage='''peaclock -i <path/to/reads> [options]
        peaclock -c <config.yaml>''')

    io_group = parser.add_argument_group('input output options')
    io_group.add_argument('-c',"--configfile",help="Config file with PEAClock run settings",dest="configfile")
    io_group.add_argument('-i','--read-path',help="Path to the directory containing fastq files",dest="read_path")
    io_group.add_argument('-o','--output-prefix', action="store",help="Output prefix. Default: peaclock_<species>_<date>")
    io_group.add_argument('--outdir', action="store",help="Output directory. Default: current working directory")
    io_group.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    

    barcode_group = parser.add_argument_group('barcode options')
    barcode_group.add_argument('-b','--barcodes-csv',help="CSV file describing which barcodes were used on which sample",dest="barcodes_csv")
    barcode_group.add_argument('-k','--barcode-kit',help="Indicates which barcode kit was used. Default: native. Options: native, rapid, pcr, all",dest="barcode_kit")

    demux_group = parser.add_argument_group('demultiplexing options')
    demux_group.add_argument('--demultiplex',action="store_true",help="Indicates that your reads have not been demultiplexed and will run guppy demultiplex on your provided read directory",dest="demultiplex")
    demux_group.add_argument('--path-to-guppy',action="store",help="Path to guppy_barcoder executable",dest="path_to_guppy")

    run_group = parser.add_argument_group('run options')
    run_group.add_argument('-s',"--species", action="store",help="Indicate which species is being sequenced. Options: mus, apodemus", dest="species")
    run_group.add_argument("-r","--report",action="store_true",help="Generate markdown report of estimated age")
    
    misc_group = parser.add_argument_group('misc options')
    misc_group.add_argument('-t', '--threads', action='store',type=int,help="Number of threads")
    misc_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")
    misc_group.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    misc_group.add_argument("-v","--version", action='version', version=f"peaclock {__version__}")

    """
    Exit with help menu if no args supplied
    """

    args = parser.parse_args(sysargs)
    
    """
    Initialising dicts
    """

    config = qcfunk.get_defaults()

    configfile = qcfunk.look_for_config(args.configfile,cwd,config)

    # if a yaml file is detected, add everything in it to the config dict
    if configfile:
        qcfunk.parse_yaml_file(configfile, config)
    else:
        if len(sysargs)<1: 
            parser.print_help()
            sys.exit(0)
    
    """
    Get outdir, tempdir and the data
    """
    # default output dir
    qcfunk.get_outdir(args.outdir,args.output_prefix,cwd,config)

    # specifying temp directory, outdir if no_temp (tempdir becomes working dir)
    tempdir = qcfunk.get_temp_dir(args.tempdir, args.no_temp,cwd,config)

    # get data for a particular species, and get species
    qcfunk.get_package_data(thisdir, args.species, config)

    config["cpg_header"] = qcfunk.make_cpg_header(config["cpg_sites"])

    # add min and max read lengths to the config
    qcfunk.get_read_length_filter(config)

    # looks for basecalled directory
    qcfunk.look_for_basecalled_reads(args.read_path,cwd,config)
    
    # looks for the csv file saying which barcodes in sample
    qcfunk.look_for_barcodes_csv(args.barcodes_csv,cwd,config)

    """
    Configure whether guppy barcoder needs to be run
    """

    qcfunk.look_for_guppy_barcoder(args.demultiplex,args.path_to_guppy,cwd,config)


    # don't run in quiet mode if verbose specified
    if args.verbose:
        quiet_mode = False
        config["log_string"] = ""
    else:
        quiet_mode = True
        lh_path = os.path.realpath(lh.__file__)
        config["log_string"] = f"--quiet --log-handler-script {lh_path} "

    qcfunk.add_arg_to_config("threads",args.threads,config)
    
    try:
        config["threads"]= int(config["threads"])
    except:
        sys.stderr.write(qcfunk.cyan('Error: Please specifiy an integer for variable `threads`.\n'))
        sys.exit(-1)
    threads = config["threads"]

    print(f"Number of threads: {threads}\n")

    # find the master Snakefile
    snakefile = qcfunk.get_snakefile(thisdir)

    if args.verbose:
        print("\n**** CONFIG ****")
        for k in sorted(config):
            print(qcfunk.green(k), config[k])

        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                        workdir=tempdir,config=config, cores=threads,lock=False
                                        )
    else:
        logger = custom_logger.Logger()
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=tempdir,
                                    config=config, cores=threads,lock=False,quiet=True,log_handler=logger.log_handler
                                    )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1

if __name__ == '__main__':
    main()