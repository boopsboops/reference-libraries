#!/usr/bin/env Rscript

# R script to make a 12S (Miya fragment) reference database for UK dishes
# downloads all mtDNA sequence data from GenBank, for a provided list of species 
# then, it uses a hidden markov model to pull out the Miya fragment from all that mtDNA data
# then, it queries NCBI for those accessions, and retrieves full metadata for them to allow better curation of reference database

# Rupert A. Collins :: 06/08/18

# load packages
# until changes come through to CRAN, need to run dev version of traits and install with sudo R in the terminal 
# either 'install_github("ropensci/traits")' or for root lib locations:
# sudo R
# library(devtools); withr::with_libpaths(new = "/usr/local/lib/R/site-library/", install_github("ropensci/traits"))

library("tidyverse")
library("rentrez")
library("rfishbase")
library("traits")
library("parallel")
library("ape")

# function for removing duplicate haplotypes
# source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/hapCollapse.R")
# function for running hmmer in R
source("run_hmmer.R")

### ### ### ###

# timing the code
#start_time <- Sys.time()

# load up the uk species table
uk.species.table <- read_csv(file="../species/uk_species_table.csv")

# make a query for genbank
range <- "1:20000" # includes mt genomes
query <- paste0("(", uk.species.table$sciName, "[ORGN] AND mitochondrion[ALL] AND ", range, "[SLEN]) OR (", uk.species.table$sciName, "[ORGN] AND mitochondrial[ALL] AND ", range, "[SLEN])")

# run the search for accessions with rentrez - takes about 10 mins
search.res <- mcmapply(FUN=function(x) entrez_search(db="nuccore", term=x, retmax=1000), query, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=8)
# if parallel package not available, use lapply (slower)
# search.res <- lapply(query, function(x) entrez_search(db="nuccore", term=x, retmax=1000))

# removed the zero counts
search.res.nz <- search.res[which(lapply(search.res, function(x) x$count) > 0)]

# get IDs and remove dups
search.ids <- unique(unlist(lapply(search.res.nz, function(x) x$ids)))

# chunk up into 300s to stop server from rejecting request 
chunk <- 300
id.split <- unname(split(search.ids, ceiling(seq_along(search.ids)/chunk)))

# download with ape (fast)
ncbi.all <- mcmapply(FUN=function(x) read.GenBank(x, species.names=FALSE), id.split, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=8)

# write out a temporary file (it's about 70 MB) for hmmer to pull out the 12S Miya sequences 
lapply(ncbi.all, write.dna, file="../temp/ncbi_uk.fas", format="fasta", append=TRUE, colw=99999)

# now run hmmer
# need to have "hmmer" and "biosquid" installed 
# if not, run 'sudo apt install hmmer biosquid'
# assumes the hidden markov model is located in hmms/ directory and is named '$prefix.hmm'
# returns a DNAbin object of the sequences matched by hmmer 
prefix <- "12s.miya.primers"
prefix <- "12s.miya.noprimers"
dat.frag <- run_hmmer(dir="../temp", infile="ncbi_uk.fas", prefix=prefix, evalue="4e-10")

# delete that temp fasta file
#file.remove("../temp/ncbi_uk.fas")

# now for the same sequences, get the tabular data from NCBI using 'ncbi_byid' to make a proper reference database
chunk <- 100
chunk.frag <- unname(split(names(dat.frag), ceiling(seq_along(names(dat.frag))/chunk)))
ncbi.frag <- mcmapply(FUN=ncbi_byid, chunk.frag, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=8)
# if parallel package not available, use lapply (slower)
# ncbi.frag <- lapply(chunk.frag, ncbi_byid)

# join all the data sets
frag.df <- as.tibble(bind_rows(ncbi.frag))

# timing the code
#end_time <- Sys.time()
#end_time - start_time

# write out to check
write_csv(frag.df, path="../references/uk-fishes-miya-12s-noprimers.csv")

# make into a fasta file if needed 
#dtmp <- strsplit(ncbi.df$sequence, "")
#names(dtmp) <- ncbi.df$acc_no
#dat <- as.DNAbin(dtmp)
#print(dat)
#write.dna(dat, file="../temp/ncbi_uk.fas", format="fasta", colw=99999)
