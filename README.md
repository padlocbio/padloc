<a href="https://github.com/padlocbio/padloc/LICENSE" alt="License"><img src="https://img.shields.io/github/license/padlocbio/padloc" /></a> <a href="https://github.com/padlocbio/padloc/" alt="Last update"><img src="https://img.shields.io/github/last-commit/padlocbio/padloc?label=last%20update" /></a> <a href="https://github.com/padlocbio/padloc/releases" alt="Release"><img src="https://img.shields.io/github/v/release/padlocbio/padloc" /></a> 

# PADLOC: Prokaryotic Antiviral Defence LOCator

## About

PADLOC is a software tool for identifying antiviral defence systems in prokaryotic genomes. PADLOC screens genomes against a database of HMMs and system classifications to find and annotate defence systems based on sequence homology and genetic architecture. PADLOC can be installed and used via the command line or via our [web app](https://padloc.otago.ac.nz).

## Citation

Manuscript in preparation

## Installation

PADLOC can be installed by cloning or downloading this github repository and running the setup script.

```bash
# Clone repo to $HOME
git clone https://github.com/padlocbio/padloc $HOME/padloc
# Add to $PATH
export PATH="$HOME/padloc:$PATH"
# Run setup
padloc --setup
```

The setup script installs software dependencies via [Homebrew](https://brew.sh/) and R packages from [CRAN](https://cran.r-project.org/) if they are not already installed. It also downloads and compiles the default database used by PADLOC from [PADLOC-DB](https://github.com/padlocbio/padloc-db).

## Examples

```bash
# BASIC: Search an amino acid fasta file with accompanying GFF annotations
padloc --faa genome.faa --gff features.faa
```

```bash
# BASIC: Search a nucleic acid fasta file, identifying CDS with prodigal
padloc --fna genome.fna
```

```bash
# INTERMEDIATE: Use multiple cpus and save output to a different directory
padloc --faa genome.faa --gff features.faa --outdir path_to_output --cpu 4
```

```bash
# ADVANCED: Use your own HMMs and system models
padloc --faa genome.faa --gff genome.gff --data path_to_data
```

## Options

```
General:
    --help        Print this help message
    --version     Print version information
    --citation    Print citation information
    --debug       Run with debug messages
Setup:
    --setup       Set up a fresh padloc install
    --bootstrap   Install dependencies (called by --setup)
    --updatedb    Install latest database (called by --setup)
Input:
    --faa [f]     Amino acid FASTA file
    --gff [f]     GFF file (only valid with [--faa])
    --fna [f]     Nucleic acid FASTA file
Output:
    --outdir [d]  Output directory
Optional:
    --data [d]    Data directory
    --cpu [n]     Use [n] CPUs (default '1')
    --raw-out     Include a summarised raw output file
```

## Output

| File           | Description                                         |
| -------------- | --------------------------------------------------- |
| *.domtblout    | Domain table file generated by HMMER.               |
| *_prodigal.faa | Amino acid FASTA file generated by prodigal.        |
| *_prodigal.gff | GFF annotation file generated by prodigal.          |
| *_padloc.csv   | PADLOC output file for identified defence systems.  |
| *_padloc.gff   | GFF annotation file for identified defence systems. |

## PADLOC-DB

The HMMs and defence system models used by PADLOC are available from the repository [PADLOC-DB](https://github.com/leightonpayne/padloc-db). This data is downloaded and compiled automatically when `padloc --setup` or `padloc --updatedb` is run. Alternatively, a custom database can be specified with `--data [d]`, refer to [PADLOC-DB](https://github.com/leightonpayne/padloc-db) for more information about configuring a custom database.

## FAQ

- **What are the requirements for an FAA/GFF file pair as input?**

  The GFF file should conform to the [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md). Each sequence in the FAA file is matched to an entry in the GFF file based on its ID attribute e.g. for the following sequence:

  ```
  >WP_000000001.1 molybdopterin-dependent oxidoreductase, partial [Escherichia coli]
  AAAAAAAGLSVPGVARAVLVSRKPSNGIKAPCRFCGTGCGVLVGTQQGRVVACQGDPDAPVNRGLNCIKG
  YFLPKIMYGKDRLTQPLLRMKNGKYDKEGEFTPITWDQAFDVMEEKFKTALKEKGPESIGMFGSGQWTIW
  EGYAASKLFKAGFRSNNIDPNARHCMASAVVGFMRTFGMDEPMGCYDDIEQADAFVLWGANMAEMHPILW
  SRITNRRLSN
  ```

  The corresponding entry in the GFF file should contain an ID attribute of the form:

  `ID=WP_000000001.1 ` *or* ``ID=cds-WP_000000001.1 ``

  FAA/GFF combinations that are known to work 'out-of-the-box' are from genomes annotated with:

  - [NCBI's prokaryotic genome annotation pipeline](https://doi.org/10.1093/nar/gkw569) (i.e. genomes from [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) and [GenBank](https://www.ncbi.nlm.nih.gov/genbank/)) 
  - [JGI's IMG annotation pipeline](https://img.jgi.doe.gov/docs/pipelineV5/) (i.e. genomes from [IMG](https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=TreeFile&page=domain&domain=all))
  - [Prokka](https://github.com/tseemann/prokka)

- **Why are there parsing failures when using a GFF file from prokka?**

  The following warning may be thrown when using a GFF file generated by prokka:

  ```
  Warning: 46324 parsing failures.
   row col  expected    actual         file
  2612  -- 9 columns 1 columns 'prokka.gff'
  2613  -- 9 columns 1 columns 'prokka.gff'
  2614  -- 9 columns 1 columns 'prokka.gff'
  2615  -- 9 columns 1 columns 'prokka.gff'
  2616  -- 9 columns 1 columns 'prokka.gff'
  .... ... ......... ......... ............
  See problems(...) for more details.
  ```

  This is because these GFF files are appended with the contig sequences of the annotated genome. This warning can be avoided by removing the contig sequences from the GFF file with:

  ```bash
  sed '/^##FASTA/Q' prokka.gff > nosequence.gff
  ```

## Issues

Bugs and feature requests can be submitted to the [Issues tab](https://github.com/leightonpayne/padloc/issues).

## Dependencies

These dependencies are installed automatically during setup.

### Mandatory

- **R**  
  *R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.*
  - **tidyverse**  
    *Wickham et al. (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686.*
  - **yaml**  
    *Stephens, J., et al. (2020). yaml: Methods to Convert R Data to YAML and Back. https://CRAN.R-project.org/package=yaml.*
  - **getopt**  
    *Davis, T., et al. (2019). getopt: C-Like 'getopt' Behavior. https://CRAN.R-project.org/package=getopt.*
- **HMMER3**  
*Finn, R.D., Clements, J., and Eddy, S.R. (2011). HMMER web server: interactive sequence similarity searching. Nucleic Acids Res 39, W29–W37.*

### Optional

- **Prodigal**  
  *Hyatt, D., Chen, GL., Locascio, P.F., Land, M.L., Larimer, F.W., and Hauser, L.J. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119.*

## License

This software and data is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
