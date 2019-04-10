#!/usr/bin/env Rscript

# R script to make reference databases for UK fishes for multiple markers
# then, it uses a hidden markov model to pull out the fragment of interest from all that mtDNA data
# then, it queries NCBI/BOLD for those accessions, and retrieves full metadata for them to allow better curation of reference database
# output is a csv dataframe of all accessions with data for all primer sets if present

# load functions and libs
source("funs.R")


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
#"coi.lerayxt.primers",
"coi.lerayxt.noprimers",
#"coi.seamid.primers",
"coi.seamid.noprimers",
#"coi.seashort.primers",
"coi.seashort.noprimers",
#"coi.ward.primers",
"coi.ward.noprimers",
#"12s.miya.primers",
"12s.miya.noprimers",
#"12s.riaz.primers",
"12s.riaz.noprimers",
#"12s.valentini.primers",
"12s.valentini.noprimers",
#"12s.taberlet.primers",
"12s.taberlet.noprimers",
"16s.berry.noprimers",
"cytb.minamoto.noprimers")


# run hmmer (takes about 5 mins)
dat.frag.all <- lapply(prefixes.all, function(x) run_hmmer3(dir="../temp", infile="mtdna-uk.fas", prefix=x, evalue="1", coords="env"))

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
ncbi.frag <- mcmapply(FUN=ncbi_byid, chunk.frag, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=2)# mc.cores=1 is the safest option, but try extra cores to speed up if there are no errors

# check for errors (should all be "data.frame")
table(sapply(ncbi.frag,class))


##
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
dat.frag.flat <- lapply(dat.frag.all, function(x) mcmapply(str_flatten, as.character(x), mc.cores=8, SIMPLIFY=TRUE,USE.NAMES=TRUE))

# turn each into a dataframe
dat.frag.df <- lapply(dat.frag.flat, function(x) tibble(names=names(x), seqs=unlist(x), lengthFrag=str_length(seqs)))

# rename each df with names of the fragment
dat.frag.df <- mapply(function(x,y,z) dplyr::rename(x,dbid=names, !!y:=seqs, !!z:=lengthFrag), dat.frag.df, paste("nucleotidesFrag",names(dat.frag.df),sep="."), paste("lengthFrag",names(dat.frag.df),sep="."), SIMPLIFY=FALSE)

# merge all the data frames 
dat.frag.merged <- dat.frag.df %>% purrr::reduce(full_join, by="dbid")

# join with the metadata dataframe
dbs.merged.all <- dplyr::left_join(dbs.merged.all,dat.frag.merged,by="dbid")


## 
# get the proper fishbase taxonomy, not the ncbi nonsense
data(fishbase)

# make a binomial scientific name
dbs.merged.all %<>% mutate(sciNameBinomen=apply(str_split_fixed(sciNameOrig, " ", 3)[,1:2], 1, paste, collapse=" "))

# fix by hand all the binomials that aren't in fishbase (see below)
dbs.merged.all$sciNameBinomen[which(dbs.merged.all$sciNameBinomen=="Xenocypris argentea")] <- "Xenocypris macrolepis"
dbs.merged.all$sciNameBinomen[which(dbs.merged.all$sciNameBinomen=="Gobio balcanicus")] <- "Gobio gobio"
dbs.merged.all$sciNameBinomen[which(dbs.merged.all$sciNameBinomen=="Sebastes marinus")] <- "Sebastes norvegicus"
dbs.merged.all$sciNameBinomen[which(dbs.merged.all$sciNameBinomen=="Chelon ramada")] <- "Liza ramada"
dbs.merged.all$sciNameBinomen[which(dbs.merged.all$sciNameBinomen=="Hippocampus ramulosus")] <- "Hippocampus guttulatus"
dbs.merged.all$sciNameBinomen[which(dbs.merged.all$sciNameBinomen=="Cichlasoma portalagrense")] <- "Cichlasoma portalegrense"
dbs.merged.all$sciNameBinomen[which(dbs.merged.all$sciNameBinomen=="Hirundichthys volador")] <- "Exocoetus volador"


# get unique species names from db output
u.sci <- unique(dbs.merged.all$sciNameBinomen)
# validate using fishbase and select the first match
v.sci <- mclapply(u.sci, validate_names, server="fishbase", mc.cores=8)
v.sci <- lapply(v.sci, function(x) x[1])

# potential problem - search for na responses not in fishbase synonyms, and fix by hand (see above)
u.sci[which(lapply(v.sci, is.na)==TRUE)]

##
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

# clean up and remove trinomials AGAIN
dbs.merged %<>% mutate(nucleotides=str_to_lower(nucleotides)) %>% 
    mutate(sciNameValid=apply(str_split_fixed(sciNameValid, " ", 3)[,1:2], 1, paste, collapse=" ")) 

# find names that are in dbs.merged, but not in the uk species table (valid names only)
extras <- unique(dbs.merged$sciNameValid)[!unique(dbs.merged$sciNameValid) %in% uk.species.table$sciName[uk.species.table$synonym==FALSE]]
print(sort(extras))

# if theses are okay, drop them from the table
dbs.merged %<>% filter(!sciNameValid %in% extras)

# check the FB synonyms! Are these okay?
unique(paste(dbs.merged$sciNameValid[which(dbs.merged$sciNameValid != dbs.merged$sciNameOrig)], dbs.merged$sciNameOrig[which(dbs.merged$sciNameValid != dbs.merged$sciNameOrig)], sep=" | "))


##
# now clean up to make the table human readable
# drop the DNA fragments and reorder the columns 
dbs.merged.info <- dbs.merged %>% select(-matches("Frag")) %>% 
   select(source,dbid,gbAccession,sciNameValid,subphylum,class,order,family,genus,sciNameBinomen,sciNameOrig,speciesCodeFishbase,
    country,catalogNumber,institutionCode,decimalLatitude,decimalLongitude,publishedAs,publishedIn,publishedBy,
    date,notesGenBank,length,nucleotides)

# make a data frame of just the sequence data
dbs.merged.seqs <- dbs.merged %>% select(matches("Frag|dbid"))

# remerge with the reorganised dataframe and remove any ids with no nucleotides
dbs.merged.final <- left_join(dbs.merged.info,dbs.merged.seqs,by="dbid") %>%
    arrange(class,order,family,genus,sciNameValid) %>% 
    filter(!is.na(nucleotides))

# take a look 
glimpse(dbs.merged.final)
names(dbs.merged.final)

##
# write out a gzipped file (orig is too big for github)
write_csv(dbs.merged.final, path=gzfile("../references/uk-fish-references.csv.gz"), na="")
write_csv(dbs.merged.final, path="../temp/uk-fish-references.csv", na="")


# to write out a fasta
#filter(dbs.merged.final, !is.na(nucleotidesFrag.coi.lerayxt.noprimers))
#write.FASTA(tab2fas(df=filter(dbs.merged.final, !is.na(nucleotidesFrag.coi.lerayxt.noprimers)), seqcol="nucleotidesFrag.coi.lerayxt.noprimers", namecol="dbid"), file=paste0("../temp/uk-fishes.",prefix,".fas"))
