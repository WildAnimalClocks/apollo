
#!/usr/bin/env python3

import os
import argparse
import csv 
import sys
from Bio import SeqIO
from datetime import datetime 
from datetime import date
import tempfile
import pkg_resources
import yaml
import subprocess

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m'
CYAN = '\u001b[36m'
DIM = '\033[2m'

def get_defaults():
    default_dict = {
        "no_temp":False,
        "demultiplex":False,
        "path_to_guppy":False,
        "barcode_kit":"native",
        "output_prefix":"apollo",
        "species":"apodemus",
        "barcodes":"native",
        "barcodes_csv":"",
        "configfile":False,
        "allowed_species":["apodemus","mus","phalacrocorax"],
        "force":True,
        "threads":1
        }
    return default_dict

def add_arg_to_config(key,arg,config):
    if arg:
        config[key] = arg

def look_for_config(configfile_arg,cwd, config):
    configfile = ""
    if configfile_arg:
        expanded_path = os.path.expanduser(configfile_arg)
        configfile = os.path.join(cwd,expanded_path)

        path_to_file = os.path.abspath(os.path.dirname(configfile))
        config["path_to_config"] = path_to_file

        if os.path.isfile(configfile):
            config["configfile"] = configfile
        else:
            sys.stderr.write(cyan(f'Error: cannot find configfile at {configfile}\n'))
            sys.exit(-1)

    elif config["configfile"]:
        expanded_path = os.path.expanduser(config["configfile"])
        configfile = os.path.join(cwd,expanded_path)

        path_to_file = os.path.abspath(os.path.dirname(configfile))
        config["path_to_config"] = path_to_file

        if os.path.isfile(configfile):
            config["configfile"] = configfile
        else:
            configfile = ""
            config["configfile"] = ""
            print(green(f'Note: cannot find configfile input'))
    else:
        configfile = os.path.join(cwd, "config.yaml")
        
        if os.path.isfile(configfile):
            config["configfile"] = configfile
            config["path_to_config"] = cwd
        else:
            configfile = ""
            config["configfile"] = ""

    return configfile

def parse_yaml_file(configfile,config):
    with open(configfile,"r") as f:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
        for key in input_config:
            snakecase_key = key.replace("-","_").lstrip("-")
            config[snakecase_key] = input_config[key]

def make_timestamped_outdir(cwd,outdir,config):

    output_prefix = config["output_prefix"]
    split_prefix = output_prefix.split("_")
    if split_prefix[-1].startswith("20"):
        output_prefix = '_'.join(split_prefix[:-1])
    config["output_prefix"] = output_prefix
    species = config["species"]
    timestamp = str(datetime.now().isoformat(timespec='milliseconds')).replace(":","").replace(".","").replace("T","-")
    outdir = os.path.join(cwd, f"{output_prefix}_{species}_{timestamp}")
    rel_outdir = os.path.join(".",timestamp)

    return outdir, rel_outdir

def get_outdir(outdir_arg,output_prefix_arg,cwd,config):
    outdir = ''
    
    add_arg_to_config("output_prefix",output_prefix_arg, config)

    if outdir_arg:
        expanded_path = os.path.expanduser(outdir_arg)
        outdir = os.path.join(cwd,expanded_path)
        rel_outdir = os.path.relpath(outdir, cwd) 

    else:
        outdir, rel_outdir = make_timestamped_outdir(cwd,outdir,config)
    
    today = date.today()
    d = today.strftime("%Y-%m-%d")
    output_prefix = config["output_prefix"]
    split_prefix = output_prefix.split("_")
    if split_prefix[-1].startswith("20"):
        output_prefix = '_'.join(split_prefix[:-1])
    config["output_prefix"] = f"{output_prefix}_{d}"

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    print(green(f"Output dir:") + f" {outdir}")
    config["outdir"] = outdir 
    config["rel_outdir"] = os.path.join(".",rel_outdir) 

def get_temp_dir(tempdir_arg,no_temp_arg, cwd,config):
    tempdir = ''
    outdir = config["outdir"]
    if no_temp_arg:
        print(green(f"--no-temp:") + f" All intermediate files will be written to {outdir}")
        tempdir = outdir
        config["no_temp"] = no_temp_arg
    elif config["no_temp"]:
        print(green(f"--no-temp:") + f" All intermediate files will be written to {outdir}")
        tempdir = outdir
    elif tempdir_arg:
        expanded_path = os.path.expanduser(tempdir_arg)
        to_be_dir = os.path.join(cwd,expanded_path)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name

    elif "tempdir" in config:
        expanded_path = os.path.expanduser(config["tempdir"])
        to_be_dir = os.path.join(cwd,expanded_path)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name

    else:
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None)
        tempdir = temporary_directory.name
    
    config["tempdir"] = tempdir 
    return tempdir
    
def get_read_length_filter(config):
    lengths = []
    for record in SeqIO.parse(config["genes"],"fasta"):
        lengths.append(len(record))
    lengths = sorted(lengths)
    config["min_length"] = lengths[0]
    config["max_length"] = lengths[-1] + 200

def get_package_data(thisdir,species_arg, config):
    matrix_file = pkg_resources.resource_filename('apollo', "data/substitution_matrix.txt")
    config["matrix_file"] = matrix_file

    add_arg_to_config("species",species_arg, config)
    
    species = config["species"]

    if species in config["allowed_species"]:
        cpg_sites = pkg_resources.resource_filename('apollo', f"data/{species}/cpg_sites.csv")
        genes = pkg_resources.resource_filename('apollo', f"data/{species}/genes.fasta")
        primer_sequences = pkg_resources.resource_filename('apollo', f"data/{species}/primer_sequences.csv")

        config["cpg_sites"] = cpg_sites
        config["genes"] = genes
        config["primer_sequences"] = primer_sequences
    else:
        sys.stderr.write(cyan(f'Error: {species} specified not configured in apollo\nPlease select an alternative species\n'))
        sys.exit(-1)



def look_for_guppy_barcoder(demultiplex_arg,path_to_guppy_arg,cwd,config):
    add_arg_to_config("demultiplex", demultiplex_arg, config)
    add_arg_to_config("path_to_guppy", path_to_guppy_arg, config)
    
    if config["demultiplex"]:
        if config["path_to_guppy"]:
            expanded_path = os.path.expanduser(config["path_to_guppy"])
            if config["path_to_guppy"].endswith("guppy_barcoder"):
                path_to_guppy = os.path.join(cwd,expanded_path)
            else:
                path_to_guppy = os.path.join(cwd,expanded_path,"guppy_barcoder")
                config["path_to_guppy"] = path_to_guppy
            os_cmd = os.system(f"{path_to_guppy} -v")

            if os_cmd != 0:
                sys.stderr.write(cyan(f'Error: guppy_barcoder at {path_to_guppy} fails to run\n'))
                sys.exit(-1)
        
        else:
            os_cmd = os.system(f"guppy_barcoder -v")
            if os_cmd != 0:
                sys.stderr.write(cyan(f'Error: please provide the path to guppy_barcoder (`--path-to-guppy`), add guppy_barcoder to your path, or run demultiplexing in MinKNOW\n'))
                sys.exit(-1)
            else:
                config["path_to_guppy"]


def look_for_basecalled_reads(read_path_arg,cwd,config):

    if read_path_arg:
        expanded_path = os.path.expanduser(read_path_arg)
        read_path = os.path.join(cwd, expanded_path)
        if not os.path.exists(read_path):
            sys.stderr.write(cyan('Error: cannot find reads at {}\n'.format(read_path)))
            sys.exit(-1)
        else:
            config["read_path"] = read_path

    elif "read_path" in config:
        if config["read_path"]:
            expanded_path = os.path.expanduser(config["read_path"])
            read_path = os.path.join(config["path_to_config"], expanded_path)
            config["read_path"] = read_path
            fq_files = 0
            for r,d,f in os.walk(read_path):
                for fn in f:
                    filename = fn.lower()
                    if filename.endswith(".fastq") or filename.endswith(".fq"):
                        fq_files +=1

            if fq_files > 0:
                print(f"Found {fq_files} fastq files in the input directory")
            else:
                sys.stderr.write(cyan('Error: cannot find fastq files at {}\nPlease check your `--read-dir`'.format(read_path)))
                sys.exit(-1)
    else:
        sys.stderr.write(cyan('Error: `--read-dir` needed. Please input the path to the fastq read files either in the config file or via the command line.\n'))
        sys.exit(-1)

def look_for_barcodes_csv(barcodes_csv_arg,cwd,config):
    barcodes_csv = ""
    if barcodes_csv_arg:
    
        barcodes_csv = os.path.join(cwd,barcodes_csv_arg)
        if not os.path.exists(barcodes_csv):
            sys.stderr.write('Error: cannot find barcodes csv at {}\n'.format(barcodes_csv))
            sys.exit(-1)
    elif "barcodes_csv" in config:
        if config["barcodes_csv"]:
            expanded_path = os.path.expanduser(config["barcodes_csv"])
            barcodes_csv = os.path.join(config["path_to_config"], expanded_path)
    
    if barcodes_csv:
        print(f"Input barcodes csv file: {barcodes_csv}")
        barcodes = []
        with open(barcodes_csv, newline="") as f:
            reader = csv.DictReader(f)
            column_names = reader.fieldnames
            if "barcode" not in column_names:
                sys.stderr.write(f"Error: Barcode file missing header field `barcode`\n")
                sys.exit(-1)
            for row in reader: 
                if row["barcode"].startswith("NB") or row["barcode"].startswith("BC"):
                    barcodes.append(row["barcode"])
                else:
                    sys.stderr.write(f"Error: Please provide barcodes in the format `NB01` or `BC01`\n")
                    sys.exit(-1)
                    
            print(f"{len(barcodes)} barcodes read in from file")
            for i in barcodes:
                print(f"  - {i}")
            barcodes = ",".join(barcodes)
            config['barcodes_csv'] = barcodes_csv
            config["barcodes"] = barcodes
    else:
        config['barcodes_csv'] = ""
        config["barcodes"] = ""
        print(green(f"Note: No barcodes csv input"))

def check_barcode_kit():

    add_arg_to_config("barcode_kit",barcode_kit_arg, config)
    barcode_kit = config["barcode_kit"]
    if args.barcode_kit.lower() in ["native","pcr","rapid","all"]:
        config["barcode_set"] = args.barcode_kit.lower()
    else:
        sys.stderr.write(f"Error: Please enter a valid barcode kit: one of\n\t-native\n\t-pcr\n\t-rapid\n\t-all\n")
        sys.exit(-1)


def get_snakefile(thisdir):
    snakefile = os.path.join(thisdir, 'scripts','Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write(cyan(f'Error: cannot find Snakefile at {snakefile}\n Check installation\n'))
        sys.exit(-1)
    return snakefile

def make_cpg_header(cpg_csv):
    cpgs= ["sample"]
    with open(cpg_csv,"r") as f:
        cpg_file = csv.DictReader(f)
        for row in cpg_file:
            cpgs.append(row["gene"].lower()+ "_" + row["position"])
    cpg_string = ",".join(cpgs)
    return cpg_string

def colour(text, text_colour):
    bold_text = 'bold' in text_colour
    text_colour = text_colour.replace('bold', '')
    underline_text = 'underline' in text_colour
    text_colour = text_colour.replace('underline', '')
    text_colour = text_colour.replace('_', '')
    text_colour = text_colour.replace(' ', '')
    text_colour = text_colour.lower()
    if 'red' in text_colour:
        coloured_text = RED
    elif 'green' in text_colour:
        coloured_text = GREEN
    elif 'yellow' in text_colour:
        coloured_text = YELLOW
    elif 'dim' in text_colour:
        coloured_text = DIM
    elif 'cyan' in text_colour:
        coloured_text = 'cyan'
    else:
        coloured_text = ''
    if bold_text:
        coloured_text += BOLD
    if underline_text:
        coloured_text += UNDERLINE
    if not coloured_text:
        return text
    coloured_text += text + END_FORMATTING
    return coloured_text

def red(text):
    return RED + text + END_FORMATTING

def cyan(text):
    return CYAN + text + END_FORMATTING

def green(text):
    return GREEN + text + END_FORMATTING

def yellow(text):
    return YELLOW + text + END_FORMATTING

def bold_underline(text):
    return BOLD + UNDERLINE + text + END_FORMATTING



def preamble(v):
    print(green("""\n
                                           .__  .__          
                      _____  ______   ____ |  | |  |   ____  
                      \__  \ \____ \ /  _ \|  | |  |  /  _ \ 
                       / __ \|  |_> >  <_> )  |_|  |_(  <_> ) 
                      (____  /   __/ \____/|____/____/\____/ 
                           \/|__|                            
                **** Age Prediction Of Little's Little Organisms ****
""")+green(f"""
                                        {v}""")+green("""
                        ****************************************
                                                                
                            Aine O'Toole & Tom Little           
                                Edinburgh University          
\n"""))

