#!/usr/bin/env Rscript
# script to load up reference libraries and clean them up

# start timer
start_time <- Sys.time()

## load up the species info table
# species data
uk.species.table <- vroom::vroom("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/species/uk-species-table.csv",delim=",",num_threads=1,guess_max=99999,col_types=cols())
uk.species.table.orig <- uk.species.table
# remove synonyms
uk.species.table %<>% select(validName,class,order,family,genus,commonName,commonSpecies) %>% distinct()
# change taxonomy for some common species
uk.species.table %<>% filter(validName!="Pungitius laevis",validName!="Cottus perifretum",validName!="Atherina presbyter")
# filter common
uk.species.table.common <- uk.species.table %>% filter(commonSpecies==TRUE)


## load up the reference library
reflib.orig <- vroom::vroom("https://github.com/boopsboops/reference-libraries/raw/master/references/uk-fish-references.csv.gz",delim=",",num_threads=1,guess_max=99999,col_types=cols())


## clean
# load up the exclusions file to clean the data
exclusions <- vroom::vroom("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/references/exclusions.csv",delim=",",num_threads=1,guess_max=99999,col_types=cols())

# exclude bad seqs and clean
reflib.orig %<>% filter(!dbid %in% exclusions$dbid[exclusions$action=="REMOVE"])

# reassign taxonomy for some recently changed species
reflib.orig %<>% mutate(sciNameValid=str_replace_all(sciNameValid,"Pungitius laevis","Pungitius pungitius")) %>% 
    mutate(sciNameValid=str_replace_all(sciNameValid,"Cottus perifretum","Cottus gobio")) %>% 
    mutate(sciNameValid=str_replace_all(sciNameValid,"Atherina presbyter","Atherina boyeri"))

# remove unverified sequences
# remove any NA nucleotides
# remove mRNA
# remove cDNA
reflib.orig %<>% 
    filter(!is.na(nucleotides)) %>% 
    filter(!grepl("UNVERIFIED:",notesGenBank)) %>%
    filter(!grepl("similar to",notesGenBank)) %>% 
    filter(!grepl("mRNA",notesGenBank)) %>% 
    filter(!grepl("cDNA",notesGenBank)) %>% 
    filter(!grepl("transcribed",notesGenBank)) %>% 
    filter(!grepl("-like",notesGenBank)) 

# end timer
end_time <- Sys.time()

# write encouraging words
writeLines(paste(dim(reflib.orig)[1], "references loaded and filtered in", round(end_time-start_time,digits=1), "seconds ..."))
