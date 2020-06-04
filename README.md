# padloc: prokaryotic antiviral defence locator

<a href="https://github.com/leightonpayne/padloc/LICENSE" alt="Contributors"><img src="https://img.shields.io/github/license/leightonpayne/padloc" /></a> <a href="https://github.com/leightonpayne/padloc/" alt="Contributors"><img src="https://img.shields.io/github/last-commit/leightonpayne/padloc?label=last%20update" /></a> (shields don't work until repo is public)

padloc is a tool for indentifying antiviral defence systems in prokaryotic genomes. Defence genes are identified using a curated database of profile Hidden Markov Models, grouped into putative operons and screened against models that depict the typical genetic layout of many different defence systems.

- [Citation](#citation)
- [Installation](#installation)
- [Overview](#overview)
- [Usage](#usage)
- [Dependencies](#dependencies)

## Citation <a name="citation"> </a>

tbd

## Installation <a name="installation"> </a>

padloc can be installed by cloning or downloading this github repository.

1. Clone repo:

```bash
git clone https://github.com/leightonpayne/padloc
```

2. Add the padloc directory to the PATH in your profile of choice (.profile, .bashrc, .zshrc, etc.):

```bash
export PATH=$PATH:/path/to/cloned/repo/padloc
```

3. Test the command to make sure your system knows the location of padloc.sh

```bash
# If path was set correctly the following should output `/path/to/repo/padloc/padloc`
which padloc
```

## Overview <a name="overview"> </a>

```
padloc
├── data
│   ├── hmm_meta.xlsx           <- Contains metadata for the HMMs.
│   ├── sys_meta.xlsx           <- Contains metadata for the models.
│   ├── hmm
│   │   └── padlocDB.hmm        <- Hidden Markov Models of defence system proteins.
│   └── sys
│       └── *.yml               <- Models describing defence system architecture.
├── src
│   ├── check_packages.R        <- Checks that required R packages are installed.
│   └── padloc.R                <- Handles system identification.
├── LICENSE                     <- License information.
├── README.md                   <- This document.
└── padloc                      <- The main script run from the command line.
```

**IMPORTANT:** padloc ships with the `data/` directory compressed. This needs to be unpacked before use:

```bash
# Unpack then remove data.tar.gz
tar -xvzf data.tar.gz && rm data.tar.gz
```

## Usage <a name="usage"> </a>

```
Usage: 
    padloc [options] --fasta <f> --featbl <f> --out <d>

Required:
    --fasta <f>     
        Path to the input fasta file (*.faa) 
    --featbl <f>    
        Path to the input feature table file (*_feature_table.txt) 
    --out <d>       
        Path to the output directory <d>

Optional:
    --cpu <n>       
        Set number of cpus to <n>, the default is 1
    --raw-out       
        Include a summarised raw output file for each genome searched
    --append        
        Append the hmmer output to exisitng results in out_dir/domtblout/
    --cleanup
        Remove domtblout/ from the output directory if left over from 
        previous runs
    -h, --help      
        Display this help message
```

#### Input

padloc requires a genome in the form of an amino acid fasta file `*.faa`, the accompanying feature table file `*_feature_table.txt` of the particular genome, and a directory to store the output.

#### Output

The output of a typical padloc run contains the following:

```
output
├── domtblout
│   └── *.domtblout         <- Domain table files generated by HMMER.
├── *_pdlcout.csv           <- padloc output file.
├── *_around.csv            <- information for genes around the identified defence systems.
├── *_within.csv            <- information for genes within the identified defence systems.
└── *.gff                   <- GFF anotation file for identified defence systems.
```

If the same output directory is used when processing multiple genomes, `*.domtblout` files will accumulate in the `domtblout/` directory. If the `*.domtblout` file already exists for a particular genome, `hmmsearch` will not be run again for that genome unless `--append` is specified. In which case, the output of `hmmsearch` will be appended to the existing `*.domtblout` file. When processing multiple genomes in batches, this allows for the pipeline to be resumed if it is interrupted.

## Dependencies <a name="dependencies"> </a>

### R

> *R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.*

Install the R software environment from your preferred [CRAN mirror](https://cran.r-project.org/mirrors.html).

### R Packages

The following R packages also need to be installed.

#### Tidyverse

- **dplyr**; *Wickham, H., et al. (2019). dplyr: A Grammar of Data Manipulation. R package version 0.8.1. https://CRAN.R-project.org/package=dplyr*

- **plyr**; *Wickham, H. (2019). plyr: Tools for Splitting, Applying and Combining Data. R package version 1.8.5. https://CRAN.R-project.org/package=plyr*

- **readr**; *Wickham, H., et al. (2018). readr: Read Rectangular Text Data. R package version 1.3.1. https://CRAN.R-project.org/package=readr*

- **readxl**; *Wickham H., et al. (2019). readxl: Read Excel Files. R package version 1.3.0. https://CRAN.R-project.org/package=readxl*

- **rlang**; *Henry, L., et al. (2019). rlang: Functions for Base Types and Core R and 'Tidyverse' Features. R package version 0.4.0. https://CRAN.R-project.org/package=rlang*

- **stringr**; *Wickham, H. (2019). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.4.0. https://CRAN.R-project.org/package=stringr*

- **tidyr**; *Wickham, H., et al. (2019). tidyr: Easily Tidy Data with 'spread()' and 'gather()' Functions. R package version 0.8.3. https://CRAN.R-project.org/package=tidyr*

#### Other

- **yaml**; *Stephens, J., et al. (2018). yaml: Methods to Convert R Data to YAML and Back. R package version 2.2.0. https://CRAN.R-project.org/package=yaml*

- **getopt**; *Davis, T., et al. (2019). getopt: C-Like 'getopt' Behavior. R package version 1.20.3. https://CRAN.R-project.org/package=getopt*

### HMMER

> *Finn, R.D., Clements, J., and Eddy, S.R. (2011). HMMER web server: interactive sequence similarity searching. Nucleic Acids Res 39, W29–W37.*

HMMER can be installed using a precompiled package via [Homebrew](https://brew.sh/):

```bash
brew install hmmer
```

Alternatively, it can be compiled from the source code available at http://eddylab.org/sof

Refer to the [HMMER User's Guide](http://eddylab.org/software/hmmer/Userguide.pdf) for more information.
