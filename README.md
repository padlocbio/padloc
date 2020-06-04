padloc

======

prokaryotic antiviral defence locator

<a href="https://github.com/leightonpayne/padloc/LICENSE" alt="Contributors"><img src="https://img.shields.io/github/license/leightonpayne/padloc" /></a> <a href="https://github.com/leightonpayne/padloc/" alt="Contributors"><img src="https://img.shields.io/github/last-commit/leightonpayne/padloc?label=last%20update" /></a> <a href="https://github.com/leightonpayne/padlocDB/" alt="Contributors"><img src="https://img.shields.io/badge/download-padlocDB-blue" /></a> # shields don't work until repo is public


- [Citation](#citation)
- [Installation](#installation)
- [Usage](#usage)
- [Overview](#overview)
- [Dependencies](#dependencies)

## Citation <a name="citation"> </a>

tbd

## Installation <a name="installation"> </a>

padloc can be installed directly by cloning or downloading this github repository.

1. Clone repo:

```bash
git clone https://github.com/leightonpayne/padloc
```

2. Add the padloc directory to the PATH variable of your profile of choice (.bash_profile, .bashrc, etc.):

```bash
export PATH=$PATH:/path/to/cloned/repo/padloc
```

3. Test the command to make sure your system knows the location of padloc.sh

```bash
which padloc
```

If the PATH variable was set correctly, this should output:

```
/path/to/cloned/repo/padloc/padloc
```

## Usage <a name="usage"> </a>

```
padloc :: A tool for identifying phage defence systems in prokaryotic genomes

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

## Overview <a name="overview"> </a>

```
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
├── README.md                   <- This document.
└── padloc                      <- The main script run from the command line.
```

**N.B.** padloc ships with an empty  `data/` directory. These files need to be download from [padlocDB](https://github.com/leightonpayne/padlocDB) where they are maintained and updated separately:

1. Clone repo:

```bash
git clone https://github.com/leightonpayne/padlocDB
```

2. Compile database (include path to padloc to automatically move files):

```bash
./db-compile.sh "path/to/padloc" # e.g. "~/tools/padloc"
```

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
