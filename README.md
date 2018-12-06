# Reference Libraries
Metabarcoding reference libraries for UK fish species. This README outlines the contents of the repository and a brief description of the workflow involved in creating a metabarcoding reference library from scratch.

### Contents (A-Z)

* **`docs/`** - Output of the species coverage reports. Can be accessed at [boopsboops.github.io/reference-libraries/reference-library-tables.html](https://boopsboops.github.io/reference-libraries/reference-library-tables.html).
* **`hmms/`** - Hidden Markov models (HMMs) and fasta files of gene markers of interest. Currently these are various markers of COI and 12S, priming sites included and not included.
* **`mitogenomes/`** - Complete curated mitogenomes and their annotation files downloaded from the Miya's MitoFish website at [mitofish.aori.u-tokyo.ac.jp/](http://mitofish.aori.u-tokyo.ac.jp/).
* **`references/`** - Completed reference libraries in CSV format. Also contains an exclusions file containing GIs of questionable sequences.
* **`scripts/`** - R and shell scripts.
    - `activity-dates.csv` - table of dates that sequences were searched/downloaded/assembled
    - `funs.R` - helper functions, including script to run the HMMER programs via R
    - `Makefile` - makefile to generate the species coverage reports
    - `mito-build.sh` - align sequences and create HMMs 
    - `mito-download.sh` - downloads mitochondrial genomes and annotations from MitoFish
    - `mito-extract.R` - extracts the gene of interest using genome annotations
    - `qc.R` - quality control a reference library
    - `reference-library-tables.Rmd` - knitr file to prepare species coverage reports
    - `references-assemble.R` - extract and annotate reference libraries from ncbi/bold dumps 
    - `sequences-download.R` - pulls all the mitochondrial DNA from NCBI and BOLD for a list of species
    - `species-list-clean.R` - cleans and annotates a names list for UK fish species and their synonyms
    - `species-list-gbif.R` - generates a species names list for UK fish species from GBIF
* **`species/`** - Tables of species lists and tissue samples.
* **`temp/`** - Temporary file directory that is not committed to the repository, but needs to be created locally to run the scripts. Ignored by git.

### Workflow

The workflow comes in five steps: (1) assemble the mitogenome data and make the HMMs for each marker; (2) create, validate and annotate the list of UK species; (3) download all sequence data for every species and synonyms; (4) assemble reference libraries for specified markers from the downloaded sequence dump; and (5) create and update a markdown summary report.

All of these steps need not be carried out every time, depending on the goal of the update:

**Step 1** does not need to be repeated as the HMMs are generated and committed to the repository, so these scripts/data are only for reference purposes. 

**Step 2** does not need to be repeated unless there are taxonomic changes or a mistake that needs to be corrected. If a species "common" status needs to changed, the `uk-species-list.csv` and `uk-species-table.csv` can be edited directly (be sure to edit both).

**Steps 3/4**  needs to be repeated each time the reference library is refreshed with new data from GenBank. These steps should be carried out every few months.

**Step 5**  needs to be repeated after either the GenBank data is updated (i.e. when steps 3/4 are run), or after the tissue samples spreadsheet (`species/tissues.csv`) is updated.

**Step 6**  is an optional quality control step, that should be carried periodically, expecially when large numbers of sequences are added to the reference libraries.

More information is found in each individal script. Generally to identify potential errors, scripts should be run line-by-line in an R console such as RStudio rather than in batch from the terminal. All packages required are listed in `scripts/funs.R`, and are standard CRAN packages with the exception of traits, which needs to be installed via GitHub (sees script for more details). The programs HMMER, MAFFT, and RAxML need to be installed on your system. 

1. * Run `scripts/mito-download.sh` to download all fish mitochondrial genomes and annotations
   * Run `scripts/mito-extract.R` to extract the single genes (e.g. 12S) from the mitogenomes
   *  Run `scripts/mito-build.sh` to align genes and generate the HMMs

2. *  Run `scripts/species-list-gbif.R` to generate a species list from the UK
   *  Run `scripts/species-list-clean.R` to clean up and annotate the species list with synonyms/taxonomy etc

3. * Run `scripts/sequences-download.R` to get all mtDNA data from NCBI/BOLD for a species list

4. * Run `scripts/references-assemble.R` to subset (with HMMER), annotate and save reference libraries for a given primer set

5. * Type `make` in a terminal in the scripts directory (this creates a new version of  [boopsboops.github.io/reference-libraries/reference-library-tables.html](https://boopsboops.github.io/reference-libraries/reference-library-tables.html))

6. * [Optional] Run `scripts/qc.R` to quality control the reference sequences and add low quality sequences to the `species/exclusions.csv` file
