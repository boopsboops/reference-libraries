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

# timing the code
start_time <- Sys.time()

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

# set n cores 
# cores=1 is the safest option, but 8 cores is faster if there are no errors
cores <- 2

# break up into chunks
# longest query should be no larger than about 2000 chars
query.split  <- split(query, ceiling(seq_along(query)/28))

# collapse into strings of n species per string
query.cat <- unname(sapply(query.split, paste, collapse=" OR "))

# check max string length is < 2000 (ish)
writeLines(paste("There are",length(query.cat),"queries"))# num queries
writeLines(paste("Maximum query length is",max(sapply(query.cat, str_length)),"chars"))# max query length
writeLines("\nTesting ...")
tst <- query.cat[which(sapply(query.cat, str_length,USE.NAMES=FALSE)==max(sapply(query.cat, str_length,USE.NAMES=FALSE)))]#get longest
entrez_search(db="nuccore", term=tst[1], retmax=99999999, api_key=ncbi.key, use_history=FALSE)# test it works

writeLines(paste("\nNow running Rentrez on",cores,"cores ..."))
# run the search for accessions with rentrez

search.res <- mcmapply(FUN=function(x) entrez_search(db="nuccore", term=x, retmax=99999999, api_key=ncbi.key, use_history=FALSE), query.cat, SIMPLIFY=FALSE, USE.NAMES=FALSE, mc.cores=cores)
# if parallel package not available or NCBI API is throttling, use lapply (slower)
#search.res <- lapply(query.cat, function(x) entrez_search(db="nuccore", term=x, retmax=99999999, api_key=ncbi.key))

# check for errors - should be all false
table(sapply(search.res, function(x) grepl("Error",x)[1]))

# removed the zero counts
search.res.nz <- search.res[which(lapply(search.res, function(x) x$count) > 0)]

# get IDs and remove dups
search.ids <- unique(unlist(lapply(search.res.nz, function(x) x$ids)))

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

end_time <- Sys.time()
end_time-start_time 
