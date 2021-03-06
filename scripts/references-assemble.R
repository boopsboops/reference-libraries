#!/usr/bin/env Rscript

# R script to make reference databases for UK fishes for multiple markers
# then, it uses a hidden markov model to pull out the fragment of interest from all that mtDNA data
# then, it queries NCBI/BOLD for those accessions, and retrieves full metadata for them to allow better curation of reference database
# output is a csv dataframe of all accessions with data for all primer sets if present

## Load functions and libs
source("funs.R")
# set cores - mc.cores=1 is the safest option, but try extra cores to speed up if there are no errors
cores <- 2


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

prefixes.all <- c(
"coi.lerayxt.noprimers",
"coi.seamid.noprimers",
"coi.seashort.noprimers",
"coi.ward.noprimers",
"12s.miya.noprimers",
"12s.riaz.noprimers",
"12s.valentini.noprimers",
"12s.taberlet.noprimers",
"16s.berry.noprimers",
"cytb.minamoto.noprimers")

# run hmmer (takes about 5 mins)
dat.frag.all <- lapply(prefixes.all, function(x) run_hmmer3(dir="../temp", infile="mtdna-uk.fas", prefix=x, evalue="10", coords="env"))

# concatentate all
dat.frag.cat <- do.call(c,dat.frag.all)

# get unique names
dat.frag.names <- unique(labels(dat.frag.cat))

# separate the extracted sequences that are in GenBank or BOLD
in.bold <- dat.frag.names[dat.frag.names %in% bold.red$processidUniq]
in.gb <- dat.frag.names[!dat.frag.names %in% bold.red$processidUniq]

# now for the same sequences, get the tabular data from NCBI using 'ncbi_byid' to make a proper reference database
chunk <- 70
chunk.frag <- unname(split(in.gb, ceiling(seq_along(in.gb)/chunk)))
ncbi.frag <- mcmapply(FUN=ncbi_byid, chunk.frag, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=cores)

# check for errors (should all be "data.frame")
table(sapply(ncbi.frag,class))

# join all the data sets
frag.df <- as_tibble(bind_rows(ncbi.frag))

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
# run
bold.red %<>% filter(processidUniq %in% in.bold) %>% 
    filter(!genbank_accession %in% frag.df$gbAccession) %>% 
    mutate(source="BOLD",nucleotides=str_to_lower(nucleotides), length=as.character(str_length(nucleotides))) %>% 
    select(source,processidUniq,genbank_accession,species_name,lat,lon,country,institution_storing,catalognum,nucleotides,length) %>%
    rename(dbid=processidUniq,gbAccession=genbank_accession,sciNameOrig=species_name,decimalLatitude=lat,decimalLongitude=lon,institutionCode=institution_storing,catalogNumber=catalognum)

# merge gb and bold data
dbs.merged.all <- bind_rows(frag.df,bold.red)

# name each DNAbin object
names(dat.frag.all) <- prefixes.all

# extract nucleotides out of the DNAbin objects
dat.frag.flat <- lapply(dat.frag.all, function(x) mcmapply(str_flatten, as.character(x), mc.cores=cores, SIMPLIFY=TRUE,USE.NAMES=TRUE))

# turn each into a dataframe
dat.frag.df <- lapply(dat.frag.flat, function(x) tibble(names=names(x), seqs=unlist(x), lengthFrag=str_length(seqs)))

# rename each df with names of the fragment
dat.frag.df <- mapply(function(x,y,z) dplyr::rename(x,dbid=names, !!y:=seqs, !!z:=lengthFrag), dat.frag.df, paste("nucleotidesFrag",names(dat.frag.df),sep="."), paste("lengthFrag",names(dat.frag.df),sep="."), SIMPLIFY=FALSE)

# merge all the data frames 
dat.frag.merged <- dat.frag.df %>% purrr::reduce(full_join, by="dbid")

# join with the metadata dataframe
dbs.merged.all <- dplyr::left_join(dbs.merged.all,dat.frag.merged,by="dbid")


## Get the proper fishbase taxonomy, not the NCBI nonsense

# make a binomial scientific name - clean mess
dbs.merged.all %<>% mutate(sciNameBinomen=apply(str_split_fixed(sciNameOrig, " ", 3)[,1:2], 1, paste, collapse=" "))

# get up to date spp list and taxonomy
fishbase.synonyms <- rfishbase::synonyms(server="fishbase")
# clean up fishbase - just accepted names and synonyms
fishbase.synonyms.acc <- fishbase.synonyms %>% mutate(TaxonLevel=str_replace_all(TaxonLevel,"^species","Species")) %>% filter(Status=="accepted name" & TaxonLevel=="Species")
fishbase.synonyms.syn <- fishbase.synonyms %>% mutate(Status=str_replace_all(Status,"Synonym","synonym")) %>% filter(Status=="synonym")

# make ref of valid uk species 
uk.species.valid <- uk.species.table %>% select(fbSpecCode,validName,class,order,family,genus,commonName) %>% distinct()

# annotate with fishbase codes and valid species names
dbs.merged.all %<>% mutate(fbSpecCode=pull(fishbase.synonyms.acc,SpecCode)[match(sciNameBinomen,pull(fishbase.synonyms.acc,synonym))]) %>% 
    mutate(fbSpecCode=if_else(is.na(fbSpecCode),pull(fishbase.synonyms.syn,SpecCode)[match(sciNameBinomen,pull(fishbase.synonyms.syn,synonym))],fbSpecCode)) %>%
    mutate(sciNameValid=pull(uk.species.valid,validName)[match(fbSpecCode,pull(uk.species.valid,fbSpecCode))])

# drop missing taxa
missing <- dbs.merged.all %>% filter(is.na(sciNameValid)) %>% pull(sciNameOrig) %>% unique()
writeLines(paste("The following taxa could not be found in UK fishes database and have been dropped:",paste(missing,collapse=", ")))
dbs.merged.all %<>% filter(!is.na(sciNameValid))

# print all the species that had their names updated
writeLines("\nThe following taxa had their GenBank names updated using FishBase:")
dbs.merged.all %>% filter(sciNameOrig != sciNameValid) %>% select(sciNameOrig,sciNameBinomen,fbSpecCode,sciNameValid) %>% arrange(sciNameOrig) %>% distinct() %>% print(n=Inf)

# add taxonomy
dbs.merged.all %<>% mutate(subphylum="Vertebrata",
    class=pull(uk.species.valid,class)[match(fbSpecCode,pull(uk.species.valid,fbSpecCode))],
    order=pull(uk.species.valid,order)[match(fbSpecCode,pull(uk.species.valid,fbSpecCode))],
    family=pull(uk.species.valid,family)[match(fbSpecCode,pull(uk.species.valid,fbSpecCode))],
    genus=pull(uk.species.valid,genus)[match(fbSpecCode,pull(uk.species.valid,fbSpecCode))])

# clean up nucs
dbs.merged.all %<>% mutate(nucleotides=str_to_lower(nucleotides))


## Clean up to make the table human readable

# drop the DNA fragments and reorder the columns 
dbs.merged.info <- dbs.merged.all %>% select(-matches("Frag")) %>% 
   select(source,dbid,gbAccession,sciNameValid,subphylum,class,order,family,genus,sciNameBinomen,sciNameOrig,fbSpecCode,
    country,catalogNumber,institutionCode,decimalLatitude,decimalLongitude,publishedAs,publishedIn,publishedBy,
    date,notesGenBank,length,nucleotides)

# make a data frame of just the sequence data
dbs.merged.seqs <- dbs.merged.all %>% select(matches("Frag|dbid"))

# remerge with the reorganised dataframe and remove any ids with no nucleotides
dbs.merged.final <- left_join(dbs.merged.info,dbs.merged.seqs,by="dbid") %>%
    arrange(class,order,family,genus,sciNameValid) %>% 
    filter(!is.na(nucleotides))

# take a look 
glimpse(dbs.merged.final)
names(dbs.merged.final)


## Write out

# first, see what was added/removed from last time
#old <- read_csv("../references/uk-fish-references.csv.gz")
#old %>% filter(!dbid %in% pull(dbs.merged.final,dbid)) %>% filter(source=="GENBANK") %>% select(dbid,sciNameOrig,sciNameBinomen,sciNameValid) %>% print(n=Inf)
#dbs.merged.final %>% filter(!dbid %in% pull(old,dbid)) %>% select(dbid,sciNameOrig,sciNameBinomen,sciNameValid) %>% print(n=Inf)

# write out a gzipped file (orig is too big for github)
write_csv(dbs.merged.final, path=gzfile("../references/uk-fish-references.csv.gz"), na="")
write_csv(dbs.merged.final, path="../temp/uk-fish-references.csv", na="")

# to write out a fasta
#filter(dbs.merged.final, !is.na(nucleotidesFrag.coi.lerayxt.noprimers))
#write.FASTA(tab2fas(df=filter(dbs.merged.final, !is.na(nucleotidesFrag.coi.lerayxt.noprimers)), seqcol="nucleotidesFrag.coi.lerayxt.noprimers", namecol="dbid"), file=paste0("../temp/uk-fishes.",prefix,".fas"))
