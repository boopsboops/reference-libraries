#!/usr/bin/env Rscript
# script to load up reference libraries and clean them up

## load up the species info table
# species data
uk.species.table <- read_csv(file="https://raw.githubusercontent.com/boopsboops/reference-libraries/master/species/uk-species-table.csv")
# remove synonyms
uk.species.table %<>% filter(synonym==FALSE)
# change taxonomy for some common species
uk.species.table %<>% filter(sciName!="Pungitius laevis",sciName!="Cottus perifretum",sciName!="Atherina presbyter")
# filter common
uk.species.table.common <- uk.species.table %>% filter(commonSpecies==TRUE)


## load up the reference library
reflib <- read_csv("https://github.com/boopsboops/reference-libraries/raw/master/references/uk-fish-references.csv.gz", guess_max=100000)


## clean
# load up the exclusions file to clean the data
exclusions <- read_csv(file="https://raw.githubusercontent.com/boopsboops/reference-libraries/master/references/exclusions.csv")

# exclude bad seqs and clean
reflib %<>% filter(!dbid %in% exclusions$dbid[exclusions$action=="REMOVE"])

# reassign taxonomy for some recently changed species
reflib %<>% mutate(sciNameValid=str_replace_all(sciNameValid,"Pungitius laevis","Pungitius pungitius")) %>% 
    mutate(sciNameValid=str_replace_all(sciNameValid,"Cottus perifretum","Cottus gobio")) %>% 
    mutate(sciNameValid=str_replace_all(sciNameValid,"Atherina presbyter","Atherina boyeri"))

# remove unverified sequences
reflib %<>% filter(!grepl("UNVERIFIED:",notesGenBank))
