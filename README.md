# PEAClock
Predicted Epigenetic Age Clock

## Quick links

  * [Requirements](#requirements)
  * [Installation](#installation)
  * [Check the install worked](#check-the-install-worked)
  * [Updating peaclock](#updating-peaclock)
  * [Input options](#input-options)
  * [Usage](#usage)
  * [Updating peaclock](#updating-peaclock)
  * [Output options](#output-options)


```

                                            __                 __     
                ______   ____ ____    ____ |  |   ____   ____ |  | __ 
                \____ \_/ __ \__  \ _/ ___\|  |  /  _ \_/ ___\|  |/ / 
                |  |_> >  ___/ / __ \  \___|  |_(  <_> )  \___|    <  
                |   __/ \___ / ____  /\___ /____/\____/ \___  >__|_ \ 
                |__|                                                   
                     **** Predicted Epigenetic Age Clock****
                
                                        0.1
                        ****************************************
                                                                
                            Aine O'Toole & Tom Little           
                                Edinburgh University          

```

### Requirements

<strong>peaclock</strong> runs on MacOS and Linux. The conda environment recipe may not build on Windows (and is not supported) but <strong>peaclock</strong> can be run using the Windows subsystem for Linux.

1. Some version of conda, we use Miniconda3. Can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html)
2. Basecalled nanopore fastq files

### Installation

1. ``git clone https://github.com/WildANimalClocks/peaclock.git`` and ``cd peaclock``
2. ``conda env create -f environment.yml``
3. ``conda activate peaclock``
4. ``python setup.py install``

> Note: we recommend using peaclock in the conda environment specified in the ``environment.yml`` file as per the instructions above. If you can't use conda for some reason, dependency details can be found in the ``environment.yml`` file. 

### Check the install worked

Type (in the peaclock environment):

```
peaclock -v
```
and you should see the version number of peaclock printed


### Updating peaclock

> Note: Even if you have previously installed ``peaclock``, as it is being worked on intensively, we recommend you check for updates before running.

To update:

1. ``conda activate peaclock``
2. ``git pull`` \
pulls the latest changes from github
3. ``python setup.py install`` \
re-installs peaclock
4. ``conda env update -f environment.yml`` \
updates the conda environment 

### Troubleshooting update
- If you have previously installed peaclock using ``pip``, you will need to update peaclock in the same way (``pip install .``)
- Try ``pip uninstall peaclock`` and then re-install with `python setup.py install`

### Input options

<strong> -c / --config</strong>

PEAClock can accept a config file in yaml (or yml) format. This is a standard config format that describes the analysis you want to run.

You can provide any of the command line arguments via this config file, for instance pass the `species` or `read_path` in through via the configuation file.

Using this input option will allow the user to run similar reports again and again, without having to specify all arguments via the command line.

Note, if the same option is specified in the config file and as a command line argument, the command line argument will overwrite the config file option. 
```
peaclock -c config.yaml
```

If no config file is specified via the command line, peaclock will look for a file called `config.yaml` in the current working directory.

Example config.yaml file:

```
demultiplex: True
outdir: config_test

read_path: path/to/mus_basecalled
path_to_guppy: ~/ont-guppy-cpu/bin/guppy_barcoder

barcode_kit: native

species: mus

threads: 4
```

peaclock config file notes
- Config keys are insensitive to '-' and '_' differences. 
  Example: `--read-path` can be added to the config.yaml file as `read_path ` or `read-path`
- All command line options are available to be input in the config file 

<strong>Command line args</strong>

See the help menu of peaclock for full command line options, all are configurable in the config.yaml file too. 

### Usage

If you have a `config.yaml` file in the same directory that you’re in, all you need to write is:
```
peaclock
```
and it’ll detect the config file and run the software with all the settings given in the config file.

If the config file isn’t in the same directory as you are (or is called soemthing different) you can say:
```
peaclock -c path/to/config.yaml
```

if you don’t have a config file and just want to run the tool on the command line you can input:

```
peaclock --read-path path/to/fastq/reads \
         --demultiplex \
         --species mus \
         -t 3 \
         --path-to-guppy path/to/guppy_barcoder
```

To run demultiplexing like above, you have either have guppy installed in your path or give peaclock the path to the binary file you download from the ont community (guppy can’t be installed with the conda command because of ont rules).

If your reads are already demultiplexed (say in MinKNOW) you can input:

```
peaclock --read-path path/to/demuxed/reads \
          --species mus \
          -t 3
```

Full usage:

```
usage: peaclock -i <path/to/reads> [options]
       peaclock -c <config.yaml>

optional arguments:
  -h, --help            show this help message and exit

input output options:
  -c CONFIGFILE, --configfile CONFIGFILE
                        Config file with PEAClock run settings
  -i READ_PATH, --read-path READ_PATH
                        Path to the directory containing fastq files
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Output prefix. Default: peaclock_<species>_<date>
  --outdir OUTDIR       Output directory. Default: current working directory
  --tempdir TEMPDIR     Specify where you want the temp stuff to go. Default:
                        $TMPDIR

barcode options:
  -b BARCODES_CSV, --barcodes-csv BARCODES_CSV
                        CSV file describing which barcodes were used on which
                        sample
  -k BARCODE_KIT, --barcode-kit BARCODE_KIT
                        Indicates which barcode kit was used. Default: native.
                        Options: native, rapid, pcr, all

demultiplexing options:
  --demultiplex         Indicates that your reads have not been demultiplexed
                        and will run guppy demultiplex on your provided read
                        directory
  --path-to-guppy PATH_TO_GUPPY
                        Path to guppy_barcoder executable

run options:
  -s SPECIES, --species SPECIES
                        Indicate which species is being sequenced. Options:
                        mus, apodemus
  -r, --report          Generate markdown report of estimated age

misc options:
  -t THREADS, --threads THREADS
                        Number of threads
  --no-temp             Output all intermediate files, for dev purposes.
  --verbose             Print lots of stuff to screen
  -v, --version         show program's version number and exit

```

## Output options

Description of output peaclock directory

<strong>-o / --output-prefix and --outdir</strong>

An output prefix can be specified with `-o /--output-prefix`. If no output directory is specified, by default the directory will be a timstamped directory beginning with the `output_prefix` and the report files will be called `output_prefix`.md  (Default: peaclock) The output directory can be specified with `--outdir`, which overwrites the `--output_prefix` directory name. See example in the figure below:
