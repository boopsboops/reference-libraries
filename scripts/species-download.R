#!/usr/bin/env Rscript

# R script to make reference databases for UK fishes for multiple markers
# downloads all mtDNA sequence data from GenBank, for a provided list of species 
# then, it uses a hidden markov model to pull out the fragment of interest from all that mtDNA data
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
source("run-hmmer.R")

### ### ### ###

# timing the code (takes about 10 mins)
#start_time <- Sys.time()

# load up the uk species table
uk.species.table <- read_csv(file="../species/uk-species-table.csv")

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

# write out a temporary file (it's about 86 MB) 
file.create("../temp/mtdna-uk.fas")
lapply(ncbi.all, write.FASTA, file="../temp/mtdna-uk.fas", append=TRUE)




# times
#end_time <- Sys.time()
#end_time-start_time 

# now run hmmer
# need to have "hmmer" and "biosquid" installed 
# if not, run 'sudo apt install hmmer biosquid'
# assumes the hidden markov model is located in hmms/ directory and is named '$prefix.hmm'
# returns a DNAbin object of the sequences matched by hmmer 
prefix <- "12s.miya.primers"
prefix <- "12s.miya.noprimers"
prefix <- "coi.lerayxt.primers"
prefix <- "coi.lerayxt.noprimers"
prefix <- "coi.seamid.primers"
prefix <- "coi.seamid.noprimers"
prefix <- "coi.seashort.primers"
prefix <- "coi.seashort.noprimers"
prefix <- "12s.riaz.primers"
prefix <- "12s.riaz.noprimers"
prefix <- "12s.valentini.primers"
prefix <- "12s.valentini.noprimers"

# run hmmer
dat.frag <- run_hmmer(dir="../temp", infile="mtdna-uk.fas", prefix=prefix, evalue="4e-10")

in.bold <- labels(dat.frag)[labels(dat.frag) %in% bold.red$processid]
in.gb <- labels(dat.frag)[!labels(dat.frag) %in% bold.red$processid]



# now for the same sequences, get the tabular data from NCBI using 'ncbi_byid' to make a proper reference database
chunk <- 100
chunk.frag <- unname(split(in.gb, ceiling(seq_along(in.gb)/chunk)))
ncbi.frag <- mcmapply(FUN=ncbi_byid, chunk.frag, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=8)
# if parallel package not available, use lapply (slower)
# ncbi.frag <- lapply(chunk.frag, ncbi_byid)

# join all the data sets
frag.df <- as.tibble(bind_rows(ncbi.frag))

# remove ncbi genome and other duplicates
frag.df <- frag.df %>% filter(gi_no!="NCBI_GENOMES") %>% 
    distinct(gi_no, .keep_all=TRUE) %>% 
    mutate(acc_no=str_replace_all(acc_no,"\\.[0-9]",""), source="GENBANK") %>%
    # fix the lat_lon into decimal
    mutate(lat=paste(str_split_fixed(lat_lon, " ", 4)[,1], str_split_fixed(lat_lon, " ", 4)[,2]), lon=paste(str_split_fixed(lat_lon, " ", 4)[,3], str_split_fixed(lat_lon, " ", 4)[,4])) %>%
    mutate(lat=if_else(grepl(" N",lat), true=str_replace_all(lat," N",""), false=if_else(grepl(" S",lat), true=paste0("-",str_replace_all(lat," S","")), false=lat))) %>%
    mutate(lon=if_else(grepl(" E",lon), true=str_replace_all(lon," E",""), false=if_else(grepl(" W",lon), true=paste0("-",str_replace_all(lon," W","")), false=lon))) %>% 
    mutate(lat=str_replace_all(lat,"^ ", NA_character_), lon=str_replace_all(lon,"^ ", NA_character_)) %>%
    mutate(lat=as.numeric(lat), lon=as.numeric(lon)) %>% 
    # tidy up
    select(-taxonomy,-organelle,-keyword,-lat_lon) %>% 
    rename(sciName=taxon,notesGenBank=gene_desc,dbid=gi_no,gbAccession=acc_no,catalogNumber=specimen_voucher,
    publishedAs=paper_title,publishedIn=journal,publishedBy=first_author,date=uploaded_date,decimalLatitude=lat,decimalLongitude=lon,nucleotides=sequence)

#bold.red %>% distinct(markercode)
bold.red <- bold.red %>% filter(processid %in% in.bold) %>% 
    filter(markercode=="COI-5P") %>% # CHANGE MARKERCODE FOR 12S
    filter(!genbank_accession %in% frag.df$gbAccession) %>% 
    mutate(source="BOLD") %>% 
    select(source,processid,genbank_accession,species_name,lat,lon,country,institution_storing,catalognum,nucleotides) %>%
    rename(dbid=processid,gbAccession=genbank_accession,sciName=species_name,decimalLatitude=lat,decimalLongitude=lon,institutionCode=institution_storing,catalogNumber=catalognum)

    
# merge gb and bold
dbs.merged <- bind_rows(frag.df,bold.red)
#write_csv(dbs.merged, path="../temp/genbank_bold_merged.csv")
    

# extract DNA from the hmmer output
nucs.list <- lapply(as.character(dat.frag), str_flatten)
nucs.df <- data_frame(names=names(nucs.list), seqs=unlist(nucs.list))

# add the just fragment
dbs.merged <- dbs.merged %>% mutate(nucleotidesFrag=nucs.df$seqs[match(dbs.merged$dbid, nucs.df$names)])

# calculate fragment lengths
dbs.merged %>% mutate(lengthFrag=str_length(nucleotidesFrag))

# get the proper fishbase taxonomy, not the ncbi nonsense
data(fishbase)

# get genera from list
# correct a couple of inconsistencies in the fishbase data
dbs.merged$sciName[is.na(fishbase[match(str_split_fixed(dbs.merged$sciName, " ", 3)[,1],fishbase$Genus),]$Order)]


u.sci <- unique(apply(str_split_fixed(dbs.merged$sciName, " ", 3)[,1:2], 1, paste, collapse=" "))
v.sci <- validate_names(u.sci, limit=1)
names.df <- data_frame(orig=u.sci, validated=v.sci)

dbs.merged %>% mutate(sciNameValid=match
names.df


# subset
fishbase.sub <- fishbase[match(str_split_fixed(frag.df$taxon, " ", 3)[,1],fishbase$Genus),]
# check length
dim(fishbase.sub)[1] == dim(frag.df)[1] 

# add taxonomy to frame and get lengths of frags and add to frame
frag.df <- frag.df %>% mutate(fragLength=str_length(frag), subphylum="Vertebrata", class=fishbase.sub$Class, order=fishbase.sub$Order, family=fishbase.sub$Family, genus=fishbase.sub$Genus)

# write out
write_csv(frag.df, path=paste0("../references/uk-fishes.",prefix,".csv"))
# frag.df <- read_csv(file=paste0("../references/uk-fishes.",prefix,".csv"))

