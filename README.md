# Reference Libraries
Metabarcoding reference libraries for UK fish species

### Contents (A-Z)

* **`docs/`** - Output of the species coverage reports. Can be accessed at [boopsboops.github.io/reference-libraries/reference-library-tables.html](https://boopsboops.github.io/reference-libraries/reference-library-tables.html).
* **`hmms/`** - Hidden Markov models (HMMs) and fasta files of gene markers of interest.
* **`mitogenomes/`** - Complete curated mitogenomes and their annotation files downloaded from the Miya's MitoFish website at [mitofish.aori.u-tokyo.ac.jp/](http://mitofish.aori.u-tokyo.ac.jp/).
* **`references/`** - Completed reference libraries in CSV format.
* **`scripts/`** - R and shell scripts.
    - `funs.R` - helper functions
    - `Makefile` - makefile to generate the species coverage reports
    - `mito-build.sh` - align sequences and create HMMs 
    - `mito-download.sh` - downloads mitochondrial genomes and annotations from MitoFish
    - `mito-extract.R` - extracts the gene of interest using genome annotations
    - `qc.R` - quality control the reference library
    - `reference-library-tables.Rmd` - knitr file to prepare species coverage reports
    - `references-assemble.R` - subset from ncbi/bold dumps and annotate for each reference marker 
    - `run-hmmer.R` - script to run the HMMER programs via R
    - `sequences-download.R` - pulls all the mitochondrial DNA from NCBI and BOLD for a list of species then pulls out the subset with the fragment of interest using the HMMs 
    - `species-list-clean.R` - cleans and annotates a names list for UK fish species and their synonyms
    - `species-list-gbif.R` - generates a species names list for UK fish species from GBIF
* **`species/`** - Tables of species lists and tissue samples.
* **`temp/`** - Temporary file directory that is not committed to the repository, but needs to be created locally to run the scripts. Ignored by git.

### Workflow

The workflow comes in three parts: (1) assemble the raw materials including mitogenome data and species lists; (2) clean and annotate the species lists and download the sequences; and (3) update the summary table document. All the files produced by steps 1 and 2 are committed into the git reposito

1. * Run `mito-download.sh` to download all fish mitochondrial genomes and annotations [06/08/2017]
   * Run `mito-extract.R` to extract the single genes (e.g. 12S) from the mitogenomes [06/08/2017]
   *  Run `mito-build.sh` to align genes and generate the HMMs [12/11/2017]
   *  Run `species-list-gbif.R` to generate a species list from the UK
   *  

2. * Run `species-list-clean.R` to clean up and annotate the species list with synonyms/taxonomy etc
   * Run `sequences-download.R` to get all mtDNA data from NCBI/BOLD for a species list
   * Run `references-assemble.R` to subset (with hmmer), annotate and save reference libraries for a given primer set

3. * Type `make` in a terminal in the scripts directory (this creates a new version of  [boopsboops.github.io/reference-libraries/reference_library_tables.html](https://boopsboops.github.io/reference-libraries/reference-library-tables.html))