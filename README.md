# A generalised, dynamic DNA reference library for UK fishes
![SeaDNA Logo](docs/logo.svg)

This repository hosts a comprehensive mitochondrial DNA reference library dataset for UK fish species, derived from the NCBI GenBank and Barcode of Life BOLD databases. The dataset includes freshwater and marine species, and can be used in a variety of applications from DNA barcoding of human food products using full COI barcodes, to metabarcoding of gut or environmental samples using fragments of 12S. The library will be updated with each new GenBank release.

A species coverage report for the current MiFish 12S dataset can be found at [boopsboops.github.io/reference-libraries/reference-library-tables.html](https://boopsboops.github.io/reference-libraries/reference-library-tables.html).

This README outlines the contents of the repository and a brief description of the workflow involved in creating a metabarcoding reference library from scratch, as well instructions to simply access the data immediately. Cloning or forking the repository will allow custom modifications to be made. If an error is apparent, raise a ticket in [Issues](https://github.com/boopsboops/reference-libraries/issues) or submit a pull request.

The work is part of the NERC funded [SeaDNA Project](https://twitter.com/SeaDNAproject), and should be cited using the following DOI: [10.6084/m9.figshare.7464521](https://doi.org/10.6084/m9.figshare.7464521).

### TL;DR

If you just want the to use the final reference database, it can be downloaded at [github.com/boopsboops/reference-libraries/raw/master/references/uk-fish-references.csv.gz](https://github.com/boopsboops/reference-libraries/raw/master/references/uk-fish-references.csv.gz). The file is 6.2 MB compressed (gzip), but is 106 MB uncompressed. The readr package function `read_csv()` will automatically decompress a gz file, but if the file needs to be unpacked to disk, the WinRAR or 7-Zip software on Windows can be used. The final dataset is in tabular CSV format; follow the Gist at [gist.github.com/boopsboops/a1c790064fe0a14af5226d098645ca60](https://gist.github.com/boopsboops/a1c790064fe0a14af5226d098645ca60) to extract the region you want, and then convert to fasta format if required. The currently available regions are as follows below in Table 1. Any additional mitochondrial primer set can be trivially added.

The dataset offered above is raw, but is cleaned when the `scripts/references-load.R` script is run. Particular attention should be paid to how this operates; sequences flagged as unreliable (using phylogenetic quality control) are listed in `references/exclusions.csv` and excluded, while sequences flagged by NCBI as "unverified" are also removed. Taxonomic changes are also made, with for example, *Cottus perifretum* relabelled as *Cottus cottus*, and *Pungitius laevis* relabelled as *Pungitius pungitius*.

**Table 1: Available primer sets**

Study | Official name | Nickname | Locus
----- | ----- | ----- | -----
[Miya et al. (2015)](http://dx.doi.org/10.1098/rsos.150088) | MiFish U/E | miya | 12S
[Taberlet et al. (2018)](http://dx.doi.org/10.1093/oso/9780198767220.001.0001) | Tele02 | taberlet | 12S
[Valentini et al. (2016)](http://dx.doi.org/10.1111/mec.13428) | L1848/H1913 | valentini | 12S
[Riaz et al. (2011)](http://dx.doi.org/10.1093/nar/gkr732) | 12S-V5 | riaz | 12S
[Wangensteen et al. (2018)](http://dx.doi.org/10.7717/peerj.4705) | Leray-XT | leray | COI
SeaDNA unpublished | SeaDNA-mid | seamid | COI
SeaDNA unpublished | SeaDNA-short | seashort | COI
[Ward et al. (2005)](http://dx.doi.org/10.1098/rstb.2005.1716) | FishF1/R1 | ward | COI

### Contents (A-Z)

* **`docs/`** - Output of the species coverage reports. Can be accessed at [boopsboops.github.io/reference-libraries/reference-library-tables.html](https://boopsboops.github.io/reference-libraries/reference-library-tables.html).
* **`hmms/`** - Hidden Markov models (HMMs) and fasta files of gene markers of interest. Currently these are various markers of COI and 12S, priming sites included and not included.
* **`references/`** - Completed reference libraries in CSV format. Also contains an exclusions file containing GIs of questionable sequences.
* **`scripts/`** - R and shell scripts.
    - `activity-dates.csv` - table of dates that sequences were searched/downloaded/assembled
    - `funs.R` - helper functions, including script to run the HMMER programs via R
    - `Makefile` - makefile to generate the species coverage reports
    - `mito-build.sh` - align sequences and create HMMs 
    - `mito-download.sh` - downloads mitochondrial genomes and annotations from the MitoFish website at [mitofish.aori.u-tokyo.ac.jp/](http://mitofish.aori.u-tokyo.ac.jp/)
    - `mito-extract.R` - extracts the gene of interest using genome annotations
    - `qc.R` - quality control a reference library
    - `reference-library-tables.Rmd` - knitr file to prepare species coverage reports
    - `references-assemble.R` - extract and annotate reference libraries from ncbi/bold dumps 
    - `references-load.R` - load up and filter/clean reference library and species list, and correct the taxonomic names
    - `sequences-download.R` - pulls all the mitochondrial DNA from NCBI and BOLD for a list of species
    - `species-list-clean.R` - cleans and annotates a names list for UK fish species and their synonyms
    - `species-list-gbif.R` - generates a species names list for UK fish species from GBIF
* **`species/`** - Tables of species lists and tissue samples.
* **`temp/`** - Temporary file directory that is not committed to the repository, but needs to be created locally to run the scripts. Ignored by git.

### Workflow

The workflow comes in five steps: (1) assemble the mitogenome data and make the HMMs for each marker; (2) create, validate and annotate the list of UK species; (3) download all sequence data for every species and synonyms; (4) assemble reference libraries for specified markers from the downloaded sequence dump; and (5) create and update a markdown summary report.

1. * Run `scripts/mito-download.sh` to download all fish mitochondrial genomes and annotations
   * Run `scripts/mito-extract.R` to extract the single genes (e.g. 12S) from the mitogenomes
   *  Run `scripts/mito-build.sh` to align genes and generate the HMMs

2. *  Run `scripts/species-list-gbif.R` to generate a species list from the UK
   *  Run `scripts/species-list-clean.R` to clean up and annotate the species list with synonyms/taxonomy etc

3. * Run `scripts/sequences-download.R` to get all mtDNA data from NCBI/BOLD for a species list

4. * Run `scripts/references-assemble.R` to subset (with HMMER), annotate and save reference libraries for a given primer set

5. * Type `make` in a terminal in the scripts directory (this creates a new version of  [boopsboops.github.io/reference-libraries/reference-library-tables.html](https://boopsboops.github.io/reference-libraries/reference-library-tables.html))

6. * [Optional] Run `scripts/qc.R` to quality control the reference sequences and add low quality sequences to the `species/exclusions.csv` file

All of these steps need not be carried out every time, depending on the goal of the update:

**Step 1** does not need to be repeated as the HMMs are generated and committed to the repository, so these scripts/data are only for reference purposes. 

**Step 2** does not need to be repeated unless there are taxonomic changes or a mistake that needs to be corrected. If a species "common" status needs to be changed, the `uk-species-list.csv` and `uk-species-table.csv` can be edited directly (be sure to edit both).

**Steps 3/4**  needs to be repeated each time the reference library needs to be refreshed with new data from GenBank. These steps should be carried out for each GenBank release every two months. The script `scripts/sequences-download.R` contains code to check the current GenBank version against the version that was last used to assemble the reference library.

**Step 5**  needs to be repeated after either the GenBank data is updated (i.e. when steps 3/4 are run), or after the tissue samples spreadsheet (`species/tissues.csv`) is updated.

**Step 6**  is an optional quality control step, that should be carried periodically, expecially when large numbers of sequences are added to the reference libraries.

More information is found in each individal script. Generally to identify potential errors, scripts should be run line-by-line in an R console such as RStudio rather than in batch from the terminal. All packages required are listed in `scripts/funs.R`, and are standard CRAN packages with the exception of [ropensci/traits](https://github.com/ropensci/traits), which needs to be installed via GitHub (sees script for more details). The programs HMMER, MAFFT, and RAxML need to be installed on your system. Unfortunately, these scripts are optimised for a Unix system, and I'm unable to offer any Windows support ([Windows is now able to run Ubuntu Linux ](https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0)).

