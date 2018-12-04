#!/usr/bin/env Rscript
# Rupert A. Collins

# R script to make reference databases for UK fishes for multiple markers
# downloads all mtDNA sequence data from GenBank/BOLD, for a provided list of species 
# last run 03/12/2018

# load libs
library("tidyverse")
library("magrittr")
library("rentrez")
library("parallel")
library("ape")
library("bold")

# function for making fasta files from tables
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")
# load up your personal NCBI API key to get 10 requests per sec. This needs to be generated from your account at https://www.ncbi.nlm.nih.gov/
# DO NOT PUT THIS KEY ON GITHUB
# if you don't have one, ncbi will rate-limit your access to 3 requests per sec, and errors may occur.
source("../temp/ncbi-key.R")
# load up the new `read_GenBank` function
source("../scripts/funs.R")

# timing the code
start_time <- Sys.time()

# load up the uk species table
uk.species.table <- read_csv(file="../species/uk-species-table.csv")


## Download all GenBank sequences for species in the UK species table (including synonyms) with mtDNA

# make a query for genbank
range <- "1:20000" # includes mt genomes
query <- paste0("(", uk.species.table$sciName, "[ORGN] AND mitochondrion[ALL] AND ", range, "[SLEN]) OR (", uk.species.table$sciName, "[ORGN] AND mitochondrial[ALL] AND ", range, "[SLEN])")

# run the search for accessions with rentrez
# sometimes an error occurs, just run again
# mc.cores=1 is the safest option, but try extra cores to speed up if there are no errors
search.res <- mcmapply(FUN=function(x) entrez_search(db="nuccore", term=x, retmax=99999, api_key=ncbi.key, use_history=TRUE), query, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=1)
# if parallel package not available or NCBI API is throttling, use lapply (slower)
#search.res <- lapply(query, function(x) entrez_search(db="nuccore", term=x, retmax=99999, api_key=ncbi.key))

# check for errors
table(sapply(search.res, function(x) grepl("Error",x)[1]))

# removed the zero counts
search.res.nz <- search.res[which(lapply(search.res, function(x) x$count) > 0)]

# get IDs and remove dups
search.ids <- unique(unlist(lapply(search.res.nz, function(x) x$ids)))

# chunk up into 100s to stop server from rejecting request 
chunk <- 100
id.split <- unname(split(search.ids, ceiling(seq_along(search.ids)/chunk)))

# download with modifed ape function (fast)
ncbi.all <- mcmapply(FUN=function(x) read_GenBank(x, species.names=FALSE, api.key=ncbi.key), id.split, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=1)

# check for errors (should be all DNAbin
table(sapply(ncbi.all, class))


# write out a temporary file
file.create("../temp/mtdna-uk.fas")
lapply(ncbi.all, write.FASTA, file="../temp/mtdna-uk.fas", append=TRUE)


## Now repeat the same for the BOLD database

# chunk up the BOLD requests
chunk <- 70
bold.split <- unname(split(uk.species.table$sciName, ceiling(seq_along(uk.species.table$sciName)/chunk)))

# query BOLD and retrieve a table
# sometimes an error occurs, just run again
bold.all <- mcmapply(FUN=function(x) bold_seqspec(x,format="tsv",sepfasta=FALSE,response=FALSE), bold.split, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=8)

# check for errors (should be "data.frame" or "logical", not "character")
table(sapply(bold.all, class))

# remove the non-dataframes
bold.all <- bold.all[which(sapply(bold.all, class)=="data.frame")]

# tidy it up and join it together
bold.red <- lapply(lapply(bold.all, as_tibble), function(x) mutate_all(x,as.character))
bold.red <- bind_rows(bold.red)
bold.red %<>% 
    mutate(nucleotides=str_replace_all(nucleotides,"-",""), nucleotides=str_replace_all(nucleotides,"N",""), num_bases=nchar(nucleotides)) %>% 
    filter(num_bases > 0) %>%
    filter(institution_storing!="Mined from GenBank, NCBI") %>% 
    mutate(processidUniq=paste(processid,markercode,sep="."))

# write temp copy of the bold dump
write_csv(bold.red, path="../temp/bold-dump.csv")

# create a fasta file of BOLD
bold.fas <- tab2fas(df=bold.red, seqcol="nucleotides", namecol="processidUniq")

# add it to the GenBank file already created
write.FASTA(bold.fas, file="../temp/mtdna-uk.fas", append=TRUE)

end_time <- Sys.time()
end_time-start_time 
