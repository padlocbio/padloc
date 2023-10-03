# Instructions for pre-computing ncRNA and CRISPR array data

> [!NOTE]
> These instructions assume that you have installed PADLOC via `conda` into an environment called `padloc`.

## Identify ncRNAs with Infernal

### Install Infernal

To simplify the installation of Infernal and dependencies, we provide the script [`padloc/etc/bin/install-infernal`](https://github.com/padlocbio/padloc/blob/master/etc/bin/install-infernal). This script creates a `conda` environment according to the spec [`padloc/etc/env/env-infernal.yml`](https://github.com/padlocbio/padloc/blob/master/etc/env/env-infernal.yml).

The script should be automatically added to your `$PATH` when you install PADLOC via `conda` and have the environment active:

```bash
conda activate padloc
install-infernal
```

Down the line, if you wish to remove the `infernal` environment, including Infernal and all dependencies, you can run:

```bash
conda env remove -n infernal
```

### Run Infernal

We also provide the script [`padloc/etc/bin/run-infernal`](https://github.com/padlocbio/padloc/blob/master/etc/bin/run-infernal) to simplify execution of Infernal with parameters consistent with the PADLOC web server. See `run-infernal --help` for full usage information. For example, running Infernal on a test genome:

```bash
run-infernal --input test/GCF_004358345.1.fna --output test/GCF_004358345.1_ncrna.tblout
```

### Run PADLOC

The formatted `.tblout.formatted` output from `run-infernal` can then be passed to PADLOC via the `--ncrna` option. For example, running PADLOC on the test genomes:

```bash
padloc --fna test/GCF_004358345.1.fna --ncrna test/GCF_004358345.1_ncrna.tblout.formatted
```

## Identify CRISPR arrays with CRISPRDetect

### Install CRISPRDetect

To simplify the installation of CRISPRDetect and dependencies, we provide the script [`padloc/etc/bin/install-crisprdetect`](https://github.com/padlocbio/padloc/blob/master/etc/bin/install-crisprdetect). This script creates a `conda` environment according to the spec [`padloc/etc/env/env-crisprdetect.yml`](https://github.com/padlocbio/padloc/blob/master/etc/env/env-crisprdetect.yml). This spec includes all dependencies required by CRISPRDetect, but not CRISPRDetect itself (which is not currently available via `conda`). CRISPRDetect is instead downloaded from the GitHub repo [`davidchyou/CRISPRDetect_2.4`](https://github.com/davidchyou/CRISPRDetect_2.4) into the new environment at `${CONDA_PREFIX}/envs/crisprdetect/CRISPRDetect_2.4/`. The perl module `Parallel::ForkManager` is also installed into the environment using `cpanm`. The PADLOC web server uses a patched version of CRISPRDetect `v2.4`. The install script applies these patches, which can be found at [`padloc/etc/patch/`](https://github.com/padlocbio/padloc/tree/master/patch).

The script should be automatically added to your `$PATH` when you install PADLOC via `conda` and have the environment active:

```bash
conda activate padloc
install-crisprdetect
```

Down the line, if you wish to remove the `crisprdetect` environment, including CRISPRDetect and all dependencies, you can run:

```bash
conda env remove -n crisprdetect
```

### Run CRISPRDetect

We also provide the script [`padloc/etc/bin/run-crisprdetect`](https://github.com/padlocbio/padloc/blob/master/etc/bin/run-crisprdetect) to simplify execution of CRISPRDetect with parameters consistent with the PADLOC web server. See `run-crisprdetect --help` for full usage information. For example, running CRISPRDetect on a test genome:

```bash
run-crisprdetect --input test/GCF_004358345.1.fna --output test/GCF_004358345.1_crispr
```

### Run PADLOC

The `.gff` output from CRISPRDetect can then be passed to PADLOC via the `--crispr` option. For example, running PADLOC on a test genomes:

```bash
padloc --fna test/GCF_004358345.1.fna --crispr test/GCF_004358345.1_crispr.gff
```

If you are specifically interested in high quality detection of CRISPR-Cas systems, we recommend also trying one of the many [specialised tools designed for this purpose](https://pubmed.ncbi.nlm.nih.gov/?term=10.1093%2Fnar%2Fgkab456+OR+10.1089%2Fcrispr.2020.0059+OR+10.1093%2Fnar%2Fgky425+OR+10.1109%2Ftcbb.2017.2665542+OR+10.1186%2Fs12864-016-2627-0+OR+10.1093%2Fnar%2Fgkaa1158+OR+10.1093%2Fgigascience%2Fgiaa062+OR+10.1089%2Fcrispr.2017.0022+OR+10.1371%2Fjournal.pone.0110726+OR+10.1002%2F1873-3468.13519+OR+10.1089%2Fcrispr.2021.0021+OR+10.1093%2Fbib%2Fbbac335+OR+10.7717%2Fpeerj.11887&sort=date) (not an exhaustive list).
