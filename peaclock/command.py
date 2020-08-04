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

import pkg_resources
from . import _program


thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description='peaclock: Predicted Epigenetic Age Clock', 
    usage='''peaclock -i <path/to/reads> -b <path/to/data> [options]''')

    parser.add_argument('-i','--read-dir',help="Input the path to the reads",dest="read_dir")
    parser.add_argument('-b','--barcodes-csv',help="CSV file describing which barcodes were used on which sample",dest="barcodes_csv")

    parser.add_argument('-s',"--species", action="store",help="Indicate which species is being sequenced", dest="species")

    # parser.add_argument("-r","--report",action="store_true",help="Generate markdown report of input queries and their local trees")

    parser.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    parser.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    parser.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")
    
    parser.add_argument('-n', '--dry-run', action='store_true',help="Go through the motions but don't actually run")
    parser.add_argument('-t', '--threads', action='store',type=int,help="Number of threads")
    parser.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    parser.add_argument("-v","--version", action='version', version=f"peaclock {__version__}")

    # Exit with help menu if no args supplied
    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    # find the master Snakefile

    snakefile = os.path.join(thisdir, 'scripts','Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n Check installation'.format(snakefile))
        sys.exit(-1)
    
    # find the barcodes file
    if args.read_dir:
        read_dir = os.path.join(cwd, args.read_dir)
        if not os.path.exists(read_dir):
            sys.stderr.write('Error: cannot find reads at {}\n'.format(read_dir))
            sys.exit(-1)
        else:
            fq_files = 0
            for r,d,f in os.walk(read_dir):
                for fn in f:
                    filename = fn.lower()
                    if filename.endswith(".fastq") or filename.endswith(".fq"):
                        fq_files +=1
            if fq_files > 0:
                print(f"Found {fq_files} fastq files in the input directory")
            else:
                sys.stderr.write('Error: cannot find fastq files at {}\n'.format(read_dir))
                sys.exit(-1)
    else:
        sys.stderr.write('Error: please input the path to the fastq read files.\n')
        sys.exit(-1)

        # find the query fasta
    if args.barcodes_csv:
        barcodes_csv = os.path.join(cwd, args.barcodes_csv)
        if not os.path.exists(barcodes_csv):
            sys.stderr.write('Error: cannot find barcodes csv at {}\n'.format(barcodes_csv))
            sys.exit(-1)
        else:
            print(f"Input barcodes csv file: {barcodes_csv}")
            with open(barcodes_csv, newline="") as f:
                reader = csv.DictReader(f)
                column_names = reader.fieldnames
                if "barcode" not in column_names:
                    sys.stderr.write(f"Error: Barcode file missing header field `barcode`\n")
                    sys.exit(-1)
    else:
        barcodes_csv = ""
        print(f"No barcodes csv file input, assuming only one sample.")

    if args.species:
        if args.species.lower() in ["apodemus","mus","phalacrocorax"]:
            species = args.species.lower()

            cpg_file = os.path.join("data",species,"cpg_sites.csv")
            reference_file = os.path.join("data",species,"genes.fasta")
            primers_file = os.path.join("data",species,"primer_sequences.csv")

            cpg_sites = pkg_resources.resource_filename('peaclock', cpg_file)
            reference_fasta = pkg_resources.resource_filename('peaclock', reference_file)
            primers = pkg_resources.resource_filename('peaclock', primers_file)

            if not os.path.isfile(cpg_sites) or not os.path.isfile(reference_fasta) or not os.path.isfile(primers):
                sys.stderr.write(f'Error: cannot find data files for {species}\n Check installation')
                sys.exit(-1)
        else:
            sys.stderr.write("""
Error: please indicate a valid species name.\n
Use one of the following:
    - mus\t\t(mouse)
    - apodemus\t\t(wood mouse)
    - phalacrocorax\t(shag)
""")
            sys.exit(-1)

    # default output dir
    outdir = ''
    if args.outdir:
        rel_outdir = args.outdir #for report weaving
        outdir = os.path.join(cwd, args.outdir)
        
        if not os.path.exists(outdir):
            os.mkdir(outdir)
    else:
        timestamp = str(datetime.now().isoformat(timespec='milliseconds')).replace(":","").replace(".","").replace("T","-")
        outdir = os.path.join(cwd, timestamp)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        rel_outdir = os.path.join(".",timestamp)
    
    print(f"Output files will be written to {outdir}\n")

    # specifying temp directory
    tempdir = ''
    if args.tempdir:
        to_be_dir = os.path.join(cwd, args.tempdir)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name
    else:
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None)
        tempdir = temporary_directory.name

    # if no temp, just write everything to outdir
    if args.no_temp:
        print(f"--no-temp: All intermediate files will be written to {outdir}")
        tempdir = outdir

    # how many threads to pass
    if args.threads:
        threads = args.threads
    else:
        threads = 1
    print(f"Number of threads: {threads}\n")

    # create the config dict to pass through to the snakemake file
    config = {
        "outdir":outdir,
        "tempdir":tempdir,
        "read_dir":read_dir,
        "rel_outdir":rel_outdir,
        "species":species,
        # "input_column":args.input_column,
        # "data_column":args.data_column,
        "force":"True"
        }


    if args.report:
        config["report"] = True
    else:
        config["report"] = False
    
    # config["report_template"] =  os.path.join(thisdir, 'scripts','report_template.pmd')
    # footer_fig = pkg_resources.resource_filename('peaclock', 'data/footer.png')
    # config["footer"] = footer_fig
    
    if args.threshold:
        try:
            threshold = int(args.threshold)
            config["threshold"] = args.threshold
        except:
            sys.stderr.write('Error: threshold must be an integer\n')
            sys.exit(-1)
    else:
        config["threshold"] = "1"

    # don't run in quiet mode if verbose specified
    if args.verbose:
        quiet_mode = False
        config["quiet_mode"]="False"
    else:
        quiet_mode = True
        config["quiet_mode"]="True"

    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=True,force_incomplete=True,workdir=tempdir,
                                 config=config, cores=threads,lock=False,quiet=quiet_mode
                                 )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1

if __name__ == '__main__':
    main()