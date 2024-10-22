<a href="https://github.com/padlocbio/padloc/blob/master/LICENSE" alt="License">
	<img src="https://img.shields.io/github/license/padlocbio/padloc" />
</a>

<a href="https://anaconda.org/padlocbio/padloc" alt="Conda release">
	<img src="https://img.shields.io/conda/vn/padlocbio/padloc?color=yellow&label=conda%20release">
</a>

<a href="https://github.com/padlocbio/padloc/commits/master" alt="GitHub commits since latest release">
	<img src="https://img.shields.io/github/commits-since/padlocbio/padloc/latest?sort=semver">
</a>

<a href="https://github.com/padlocbio/padloc/" alt="Last update">
	<img src="https://img.shields.io/github/last-commit/padlocbio/padloc?label=last%20update" />
</a>

<a href="https://doi.org/10/gzgh" alt="Citations">
	<img src="https://citations.njzjz.win/10.1093/nar/gkab883">
</a>


# PADLOC: Prokaryotic Antiviral Defence LOCator

> [!IMPORTANT]
> [PADLOC](https://github.com/padlocbio/padloc) `>v2.0.0` is only compatible with [PADLOC-DB](https://github.com/padlocbio/padloc-db) `>v2.0.0` and vice-versa. After you update PADLOC, make sure to update your database by running: `padloc --db-update`.

## About

[PADLOC](https://github.com/padlocbio/padloc) is a software tool for identifying antiviral defence systems in prokaryotic genomes. [PADLOC](https://github.com/padlocbio/padloc) screens genomes against a database of HMMs and system classifications to find and annotate defence systems based on sequence homology and genetic architecture.

## Citation

If you use [PADLOC](https://github.com/padlocbio/padloc) or [PADLOC-DB](https://github.com/padlocbio/padloc-db) please cite:

> Payne, L. J., Todeschini, T. C., Wu, Y., Perry, B. J., Ronson, C. W., Fineran, P. C., Nobrega, F. L., Jackson, S. A. (2021) Identification and classification of antiviral defence systems in bacteria and archaea with PADLOC reveals new system types. *Nucleic Acids Research*, **49**, 10868-10878. doi: https://doi.org/10.1093/nar/gkab883

If you use the [PADLOC web server](https://padloc.otago.ac.nz/padloc/) please additionally cite:

> Payne, L. J., Meaden S., Mestre M. R., Palmer C., Toro N., Fineran P. C. and Jackson S. A. (2022) PADLOC: a web server for the identification of antiviral defence systems in microbial genomes. *Nucleic Acids Research*, **50**, W541-W550. doi: https://doi.org/10.1093/nar/gkac400

The HMMs and system models in [PADLOC-DB](https://github.com/padlocbio/padloc-db) were built and curated using the data and conclusions from many different sources, we encourage you to also give credit to these groups by reading their work and citing them where appropriate. References to relevant literature can be found at the [PADLOC-DB](https://github.com/padlocbio/padloc-db/blob/master/system_info.md) repository.

## Installation

### Conda

It is recommended that [PADLOC](https://github.com/padlocbio/padloc) be installed via conda.

```bash
# Install PADLOC into a new conda environment
conda create -n padloc -c conda-forge -c bioconda -c padlocbio padloc=2.0.0
# Activate the environment
conda activate padloc
# Download the latest database
padloc --db-update
```

If you're having installation issues, refer to [Issue #35](https://github.com/padlocbio/padloc/issues/35).

## Examples

```bash
# BASIC: Search an amino acid fasta file with accompanying GFF annotations
padloc --faa genome.faa --gff features.gff
```

```bash
# INTERMEDIATE: Use multiple cpus and save output to a different directory
padloc --faa genome.faa --gff features.gff --outdir path_to_output --cpu 4
```

```bash
# ADVANCED: Supply ncRNA and CRISPR array data
padloc --faa genome.faa --gff features.gff --ncrna genome.ncrna --crispr genome.crispr
```

> [!NOTE]
> Refer to [`padloc/etc/README.md`](https://github.com/padlocbio/padloc/blob/master/etc/README.md) for instructions on pre-computing ncRNA and CRISPR array data.

## Test

```bash
# Try running PADLOC on the test data provided
padloc --faa padloc/test/GCF_001688665.2.faa --gff padloc/test/GCF_001688665.2.gff
padloc --fna padloc/test/GCF_004358345.1.fna
```

## Options

```
General:
  --help            Print this help message
  --version         Print version information
  --citation        Print citation information
  --check-deps      Check that dependencies are installed
  --debug           Run with debug messages
Database:
  --db-list         List all PADLOC-DB releases
  --db-install [n]  Install specific PADLOC-DB release [n]
  --db-update       Install latest PADLOC-DB release
  --db-version      Print database version information
Input:
  --faa [f]         Amino acid FASTA file (only valid with [--gff])
  --gff [f]         GFF file (only valid with [--faa])
  --fna [f]         Nucleic acid FASTA file
  --crispr [f]      CRISPRDetect output file containing array data
  --ncrna [f]       Infernal output file containing ncRNA data
Output:
  --outdir [d]      Output directory
Optional:
  --data [d]        Data directory
  --cpu [n]         Use [n] CPUs (default '1')
  --fix-prodigal    Set this flag when providing an FAA and GFF file
                    generated with prodigal to force fixing of sequence IDs
```

## Output

| Extension     | Description                                         |
| ------------- | --------------------------------------------------- |
| .domtblout    | Domain table file generated by HMMER.               |
| _prodigal.faa | Amino acid FASTA file generated by prodigal.        |
| _prodigal.gff | GFF annotation file generated by prodigal.          |
| _padloc.csv   | PADLOC output file for identified defence systems.  |
| _padloc.gff   | GFF annotation file for identified defence systems. |

## Interpreting Output

| Column               | Description                                                  |
| -------------------- | ------------------------------------------------------------ |
| `system.number`      | Distinct system number.                                      |
| `seqid`              | Sequence ID of the contig.                                   |
| `system`             | Name of the system identified.                               |
| `target.name`        | Protein ID.                                                  |
| `hmm.accession`      | PADLOC HMM accession number.                                 |
| `hmm.name`           | PADLOC HMM name.                                             |
| `protein.name`       | Defence system protein name.                                 |
| `full.seq.E.value`   | Full sequence E-value. From the HMMER Documentation: "The E-value is a measure of statistical signiﬁcance. The lower the E-value, the more signiﬁcant the hit." |
| `domain.iE.value`    | Domain E-value. From the HMMER Documentation: "If the full sequence E-value is signiﬁcant but the single best domain E-value is not, the target sequence is probably a multidomain remote homolog". This is the primary value used to filter the HMMER output for putative hits. |
| `target.coverage`    | Fraction of the target sequence aligning to the HMM.         |
| `hmm.coverage`       | Fraction of the HMM aligning to the target sequence.         |
| `start`              | Start position of the target sequence in the contig.         |
| `end`                | End position of the target sequence in the contig.           |
| `strand`             | Strand; forward (+) or reverse (-)                           |
| `target.description` | Target sequence descrition taken from the input file.        |
| `relative.position`  | Relative position of the target sequence in the contig.      |
| `contig.end`         | Relative position of the last sequence in the contig.        |
| `all.domains`        | Concatenated list of all domains identified with HMMER.      |
| `best.hits`          | Top 5 hits identified with HMMER.                            |

## PADLOC-DB

The HMMs and defence system models used by [PADLOC](https://github.com/padlocbio/padloc) are available from the [PADLOC-DB](https://github.com/leightonpayne/padloc-db) repository. The latest version of the database can be downloaded by running `padloc --db-update`. Alternatively, a custom database can be specified with `--data`, refer to [PADLOC-DB](https://github.com/leightonpayne/padloc-db) for more information about the database.

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

  `ID=WP_000000001.1` *or* ``ID=cds-WP_000000001.1 ``

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

- **Why can't I use a nucleotide FASTA file with < 100 kbp?**

  According to [Prodigal's own documentation](https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#plasmids-phages-viruses-and-other-short-sequences), sequences < 100 kbp are *"too short to gather enough statistics to predict genes well"*. To avoid issues arising from this, PADLOC won't try to run prodigal over anything < 100 kbp. 

  If you know what you're doing then you can use Prodigal or another gene prediction program to generate your own FAA and GFF files to then use with PADLOC.

- **Is the `full.seq.E.value` or `domain.iE.value` used to filter HMM hits?**

  The primary method of filtering hits is via the domain iE-value, based on the thresholds specified in [`hmm_meta.txt`](https://github.com/padlocbio/padloc-db/blob/master/hmm_meta.txt), under the column `e.val.threshold`.

- **How do I use `--ncrna` and `--crispr` to identify ncRNAs and CRISPR arrays?**

  With `--ncrna` and `--crispr`, pre-computed files from Infernal and CRISPRDetect respectively can be supplied to PADLOC to be included in the detection of Retrons and CRISPR-Cas systems. Infernal and CRISPRDetect are run automatically when using the PADLOC [web server](https://padloc.otago.ac.nz), but can also be run and supplied to the command line version.

> [!NOTE]
> Refer to [`padloc/etc/README.md`](https://github.com/padlocbio/padloc/blob/master/etc/README.md) for instructions on pre-computing ncRNA and CRISPR array data.

## Issues

Bugs and feature requests can be submitted to the [Issues tab](https://github.com/leightonpayne/padloc/issues) (see [Sample bug report](/../../issues/6)).

## Dependencies

These dependencies are automatically installed when installing PADLOC via `conda`.

- **R == 4.3.1**
  *R Core Team (2018). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.*
  - **tidyverse == 2.0.0**
    *Wickham et al. (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686.*
  - **yaml == 2.3.7**
    *Stephens, J., et al. (2020). yaml: Methods to Convert R Data to YAML and Back. https://CRAN.R-project.org/package=yaml.*
  - **getopt == 1.20.3**
    *Davis, T., et al. (2019). getopt: C-Like 'getopt' Behavior. https://CRAN.R-project.org/package=getopt.*
- **HMMER == 3.3.2**
*Finn, R.D., Clements, J., and Eddy, S.R. (2011). HMMER web server: interactive sequence similarity searching. Nucleic Acids Res 39, W29–W37.*
- **Prodigal == 2.6.3**
  *Hyatt, D., Chen, GL., Locascio, P.F., Land, M.L., Larimer, F.W., and Hauser, L.J. (2010) Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119.*

## License

This software and data is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
