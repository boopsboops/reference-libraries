#!/usr/bin/env Rscript
# Rupert A. Collins

# R script to make reference databases for UK fishes for multiple markers
# downloads all mtDNA sequence data from GenBank/BOLD, for a provided list of species 

# load functions and libs
source("funs.R")
# load up your personal NCBI API key to get 10 requests per sec. This needs to be generated from your account at https://www.ncbi.nlm.nih.gov/
# DO NOT PUT THIS KEY ON GITHUB
# if you don't have one, ncbi will rate-limit your access to 3 requests per sec, and errors may occur.
source("../temp/ncbi-key.R")

# load up the uk species table
uk.species.table <- read_csv(file="../species/uk-species-table.csv")

# check the GenBank data release number against the record of previous download
gb.version <- read.table("ftp://ftp.ncbi.nih.gov/genbank/GB_Release_Number")$V1
local.version <- read_csv(file="../references/activity-dates.csv",col_types=cols()) %>% filter(activity=="download mtDNA all uk species") %>% select(version) %>% pull(version)
writeLines(paste("\n\nGenBank is at version",gb.version))
writeLines(paste("Repository is at version",local.version,"\n"))

## Download all GenBank sequences for species in the UK species table (including synonyms) with mtDNA

# make a query for genbank
range <- "1:20000" # includes mt genomes, no bigger
gene.syns <- c("COI","12S","16S","rRNA","ribosomal","cytb","CO1","cox1","cytochrome","subunit","COB","CYB","mitochondrial","mitochondrion")
spp.list <- unique(c(pull(uk.species.table,speciesName),pull(uk.species.table,validName)))
query <- unlist(mapply(function(x) paste0("(",spp.list,"[ORGN] AND ",x,"[ALL] AND ",range,"[SLEN])"), gene.syns, SIMPLIFY=FALSE, USE.NAMES=FALSE))

# randomise the query
set.seed(42)
query <- sample(query,length(query))

# set n cores to parallel search in n threads
# cores=1 is the safest option, but more cores are faster if there are no errors
# do not try more than 10 cores (with api key)
# do not try more than 3 cores (without api key)
# important - try to run the search when server loads are lowest, i.e. at weekends or when the USA is not at work.
# should take about 1.5 h with 4 cores
cores <- 4

# break up into chunks
# longest query should be no larger than about 2500 chars - reduce chunk.size to get smaller queries
chunk.size <- 35
query.split  <- split(query, ceiling(seq_along(query)/chunk.size))

# collapse into strings of n species per string
query.cat <- unname(sapply(query.split, paste, collapse=" OR "))

# TEST and CHECK it's working
# check max string length is < 2500 (ish)
writeLines(paste("There are",length(query.cat),"queries"))# num queries
writeLines(paste("Maximum query length is",max(sapply(query.cat, str_length)),"chars"))# max query length
writeLines("\nTesting ...")
tst <- query.cat[which(sapply(query.cat, str_length,USE.NAMES=FALSE)==max(sapply(query.cat, str_length,USE.NAMES=FALSE)))]#get longest
options("scipen"=100)
entrez_search(db="nuccore", term=tst[1], retmax=100000, api_key=ncbi.key, use_history=FALSE)# test it works

# chunk queries over the n cores
queries.chunked  <- split(query.cat, ceiling(seq_along(query.cat)/cores))
length(queries.chunked)

# run NCBI search and time
writeLines(paste("\nNow running Rentrez on",cores,"cores ..."))
    start_time <- Sys.time()
search.res <- lapply(queries.chunked,entrez_search_parallel,threads=cores,key=ncbi.key)
    end_time <- Sys.time()
    end_time-start_time

# check for errors - should be all false
table(grepl("Error",search.res))

# flatten the searches
search.flat <- search.res %>% purrr::flatten()

# plot the counts of ids per search
search.flat %>% purrr::map(~{unname(.x$count)}) %>% purrr::flatten_int() %>% tibble() %>% ggplot(aes(x=.)) + geom_histogram()

# get max count
search.flat %>% purrr::map(~{unname(.x$count)}) %>% purrr::flatten_int() %>% max()

# get the ids out the the non-zero results
search.ids <- search.flat[which(search.flat %>% purrr::map(~{unname(.x$count)}) > 0)] %>% purrr::map(~{unname(.x$ids)}) %>% purrr::flatten_chr() %>% unique()

# count num unique ids
length(search.ids)


# now download ids using ape
# chunk up into 70s to stop server from rejecting request 
chunk <- 70
id.split <- unname(split(search.ids, ceiling(seq_along(search.ids)/chunk)))

# download with modifed ape function (fast)
# mc.cores=1 is the safest option, but more cores is faster if there are no errors
ncbi.all <- mcmapply(FUN=function(x) read_GenBank(x, species.names=FALSE, api.key=ncbi.key), id.split, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=cores)

# check for errors (should be all DNAbin)
table(sapply(ncbi.all, class))

# write out a temporary file
file.create("../temp/mtdna-uk.fas")
lapply(ncbi.all, write.FASTA, file="../temp/mtdna-uk.fas", append=TRUE)


## Now repeat the same for the BOLD database

# chunk up the BOLD requests
chunk <- 70
bold.split <- unname(split(spp.list, ceiling(seq_along(spp.list)/chunk)))

# query BOLD and retrieve a table
# sometimes an error occurs, just run again
bold.all <- mcmapply(FUN=function(x) bold_seqspec(x,format="tsv",sepfasta=FALSE,response=FALSE), bold.split, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=cores)

# check for errors (should be "data.frame" or "logical", not "character")
table(sapply(bold.all, class))

# remove the non-dataframes
bold.all <- bold.all[which(sapply(bold.all, class)=="data.frame")]

# tidy it up and join it together, remove duplicate records
bold.red <- lapply(lapply(bold.all, as_tibble), function(x) mutate_all(x,as.character))
bold.red <- bind_rows(bold.red)
bold.red %<>% 
    mutate(nucleotides=str_replace_all(nucleotides,"-",""), nucleotides=str_replace_all(nucleotides,"N",""), num_bases=nchar(nucleotides)) %>% 
    filter(num_bases > 0) %>%
    filter(institution_storing!="Mined from GenBank, NCBI") %>% 
    mutate(processidUniq=paste(processid,markercode,sep=".")) %>% 
    distinct(processidUniq, .keep_all=TRUE)

# write temp copy of the bold dump
write_csv(bold.red, path="../temp/bold-dump.csv")

# create a fasta file of BOLD
bold.fas <- tab2fas(df=bold.red, seqcol="nucleotides", namecol="processidUniq")

# add it to the GenBank file already created
write.FASTA(bold.fas, file="../temp/mtdna-uk.fas", append=TRUE)
