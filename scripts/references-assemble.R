#!/usr/bin/env Rscript

# R script to make reference databases for UK fishes for multiple markers
# then, it uses a hidden markov model to pull out the fragment of interest from all that mtDNA data
# then, it queries NCBI/BOLD for those accessions, and retrieves full metadata for them to allow better curation of reference database

# until changes come through to CRAN, need to run dev version of traits and install with sudo R in the terminal 
# either 'install_github("ropensci/traits")' or for root lib locations:
# sudo R
# library(devtools); withr::with_libpaths(new = "/usr/local/lib/R/site-library/", install_github("ropensci/traits"))

# load libs
library("tidyverse")
library("magrittr")
library("rfishbase")
library("traits")
library("parallel")
library("ape")

# function for running hmmer in R
source("run-hmmer.R")


## Data
# load up the uk species table
uk.species.table <- read_csv(file="../species/uk-species-table.csv")
# load the BOLD dump
bold.red <- read_csv(file="../temp/bold-dump.csv", guess_max=100000)

## Extract the frag of interest using the HMMs
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
prefix <- "12s.taberlet.primers"
prefix <- "12s.taberlet.noprimers"

# run hmmer
#dat.frag <- run_hmmer(dir="../temp", infile="mtdna-uk.fas", prefix=prefix, evalue="0.05")
dat.frag <- run_hmmer3(dir="../temp", infile="mtdna-uk.fas", prefix=prefix, evalue="0.05")

# separate the extracted sequences that are in GenBank or BOLD
in.bold <- labels(dat.frag)[labels(dat.frag) %in% bold.red$processidUniq]
in.gb <- labels(dat.frag)[!labels(dat.frag) %in% bold.red$processidUniq]

# now for the same sequences, get the tabular data from NCBI using 'ncbi_byid' to make a proper reference database
chunk <- 80
chunk.frag <- unname(split(in.gb, ceiling(seq_along(in.gb)/chunk)))
ncbi.frag <- mcmapply(FUN=ncbi_byid, chunk.frag, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=8)# sometimes this errors, so repeat a few times

# join all the data sets
frag.df <- as.tibble(bind_rows(ncbi.frag))

# from GenBank remove ncbi genome and other duplicates etc, and clean the lat/lon data
frag.df %<>% filter(gi_no!="NCBI_GENOMES") %>% 
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
    rename(sciNameOrig=taxon,notesGenBank=gene_desc,dbid=gi_no,gbAccession=acc_no,catalogNumber=specimen_voucher,
    publishedAs=paper_title,publishedIn=journal,publishedBy=first_author,date=uploaded_date,decimalLatitude=lat,decimalLongitude=lon,nucleotides=sequence)

# do the same for BOLD
#bold.red %>% distinct(markercode)
# first check if anything exists in BOLD
bold.red %>% filter(processidUniq %in% in.bold) %>% filter(!genbank_accession %in% frag.df$gbAccession)
# run
bold.red %<>% filter(processidUniq %in% in.bold) %>% 
    filter(!genbank_accession %in% frag.df$gbAccession) %>% 
    mutate(source="BOLD") %>% 
    select(source,processid,genbank_accession,species_name,lat,lon,country,institution_storing,catalognum,nucleotides) %>%
    rename(dbid=processid,gbAccession=genbank_accession,sciNameOrig=species_name,decimalLatitude=lat,decimalLongitude=lon,institutionCode=institution_storing,catalogNumber=catalognum)

# merge gb and bold data
dbs.merged.all <- bind_rows(frag.df,bold.red)

# extract DNA from the hmmer output and add the fragment to the merged db
nucs.list <- lapply(as.character(dat.frag), str_flatten)
nucs.df <- data_frame(names=names(nucs.list), seqs=unlist(nucs.list))
nucs.df %<>% mutate(names=str_split_fixed(names,"\\.",2)[,1])

dbs.merged.all %<>% mutate(nucleotidesFrag=nucs.df$seqs[match(dbs.merged.all$dbid, nucs.df$names)], lengthFrag=str_length(nucleotidesFrag))
# take a quick look at the fragment sizes
table(dbs.merged.all$lengthFrag)


## Now add the taxonomy 

# get the proper fishbase taxonomy, not the ncbi nonsense
data(fishbase)

# make a binomial scientific name
dbs.merged.all %<>% mutate(sciNameBinomen=apply(str_split_fixed(sciNameOrig, " ", 3)[,1:2], 1, paste, collapse=" "))

# get unique species names from db output
u.sci <- unique(dbs.merged.all$sciNameBinomen)
# validate using fishbase and select the first match
v.sci <- mclapply(u.sci, validate_names, server="fishbase")
v.sci <- lapply(v.sci, function(x) x[1])

# potential problem - search for na responses not in fishbase, and fix by hand :(
# go back to v.sci and repeat validate_names
u.sci[which(lapply(v.sci, is.na)==TRUE)]
u.sci[which(u.sci=="Xenocypris argentea")] <- "Xenocypris macrolepis"
u.sci[which(u.sci=="Gobio balcanicus")] <- "Gobio gobio"
u.sci[which(u.sci=="Sebastes marinus")] <- "Sebastes norvegicus"
u.sci[which(u.sci=="Chelon ramada")] <- "Liza ramada"

# replace the null values with a value from "limit=1" - not sure why these were generated
#v.sci[which(lapply(v.sci, function(x) is.null(x))==TRUE)] <- validate_names(u.sci[which(lapply(v.sci, function(x) is.null(x))==TRUE)], limit=1)

# make a df
names.df <- data_frame(orig=u.sci, validated=unlist(v.sci))

# add the validated names to the merged db
dbs.merged <- dbs.merged.all %>% mutate(sciNameValid=names.df$validated[match(dbs.merged.all$sciNameBinomen,names.df$orig)])

# create new fishbase species/genus and subset fishbase
fishbase %<>% mutate(genusSpecies=paste(Genus,Species))
fishbase.sub <-fishbase[match(dbs.merged$sciNameValid,fishbase$genusSpecies),]
# check length is okay
dim(fishbase.sub)[1] == dim(dbs.merged)[1] 

# add taxonomy 
dbs.merged %<>% mutate(subphylum="Vertebrata", class=fishbase.sub$Class, order=fishbase.sub$Order, family=fishbase.sub$Family, genus=fishbase.sub$Genus,speciesCodeFishbase=fishbase.sub$SpecCode)

# clean up and sort columns and remove trinomials AGAIN
dbs.merged %<>% mutate(nucleotides=str_to_lower(nucleotides)) %>% 
    mutate(sciNameValid=apply(str_split_fixed(sciNameValid, " ", 3)[,1:2], 1, paste, collapse=" ")) %>% 
    select(source,dbid,gbAccession,sciNameValid,subphylum,class,order,family,genus,sciNameBinomen,sciNameOrig,speciesCodeFishbase,
    country,catalogNumber,institutionCode,decimalLatitude,decimalLongitude,publishedAs,publishedIn,publishedBy,
    date,notesGenBank,length,lengthFrag,nucleotidesFrag,nucleotides) %>% 
    arrange(class,order,family,genus,sciNameValid)

# find names that are in dbs.merged, but not in the uk species table (valid names only)
extras <- unique(dbs.merged$sciNameValid)[!unique(dbs.merged$sciNameValid) %in% uk.species.table$sciName[uk.species.table$synonym==FALSE]]
print(extras)

# if theses are okay, drop them from the table
dbs.merged %<>% filter(!sciNameValid %in% extras)

# check the FB synonyms! Are these okay?
unique(paste(dbs.merged$sciNameValid[which(dbs.merged$sciNameValid != dbs.merged$sciNameOrig)], dbs.merged$sciNameOrig[which(dbs.merged$sciNameValid != dbs.merged$sciNameOrig)], sep=" | "))

# write out
#write_csv(dbs.merged, path=paste0("../references/uk-fishes.",prefix,".csv"))
write_csv(dbs.merged, path=paste0("../temp/uk-fishes.",prefix,".csv"))
# important - clear memory before re-running
rm(list=ls())

######## checking #######
oriig <- read_csv(file="../temp/libraries_29-10-18/uk-fishes.coi.seamid.noprimers.csv")

setdiff(oriig$dbid,dbs.merged$dbid)
setdiff(dbs.merged$dbid,oriig$dbid)

mis <- oriig[match(setdiff(oriig$dbid,dbs.merged$dbid), oriig$dbid),]
mis %<>% filter(!gbAccession %in% dbs.merged$gbAccession)

write_csv(mis, path="../temp/missing.csv")


bold_seqspec("Nessorhamphus ingolfianus",format="tsv",sepfasta=FALSE,response=FALSE)

which(uk.species.table$sciName=="Phoxinus phoxinus")

bold.red[which(bold.red$species_name=="Nessorhamphus ingolfianus"),]


range <- "1:20000" # includes mt genomes
query <- paste0("(", "Phoxinus phoxinus", "[ORGN] AND mitochondrion[ALL] AND ", range, "[SLEN]) OR (","Phoxinus phoxinus", "[ORGN] AND mitochondrial[ALL] AND ", range, "[SLEN])")
validate_names(species="Chelon ramada",server="fishbase")