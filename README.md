# Reference Libraries
Metabarcoding reference libraries for UK fish species

### Contents

* **`docs/`** - Output of the species coverage reports. Can be accessed at [boopsboops.github.io/reference-libraries/reference_library_tables.html](https://boopsboops.github.io/reference-libraries/reference_library_tables.html).
* **`hmms/`** - Hidden Markov models (HMMs) and fasta files of gene markers of interest.
* **`mitogenomes/`** - Complete curated mitogenomes and their annotation files downloaded from the Miya's MitoFish website at [mitofish.aori.u-tokyo.ac.jp/](http://mitofish.aori.u-tokyo.ac.jp/).
* **`references/`** - Completed reference libraries in CSV format.
* **`scripts/`** - R and shell scripts.
    - `Makefile` - makefile to generate the species coverage reports
    - `mito-build.sh` - align sequences and create HMMs 
    - `mito-download.sh` - downloads mitochondrial genomes and annotations from MitoFish
    - `mito-extract.R` - extracts the gene of interest using genome annotations
    - `reference-library-tables.Rmd` - knitr file to prepare species coverage reports
    - `run-hmmer.R` - script to run the HMMER programs via R
    - `species-download.R` - pulls all the mitochondrial DNA for a list of species then pulls out the subset with the fragment of interest using the HMMs 
    - `species-list.R` - generates an annotated species names list for UK fish species and all their synonyms  
* **`species/`** - Tables of species lists and tissue samples.
* **`temp/`** - Temporary file directory that is not committed to the repository, but needs to be created locally to run the scripts. Ignored by git.