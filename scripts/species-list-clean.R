#!/usr/bin/env Rscript

# this script cleans and annotates the uk fish species list using fishbase

# load functions and libs
source("funs.R")

# load up the combined species list (generated 20/12/2017)
uk.list <- read_csv(file="../species/uk-species-list.csv")

# get up to date spp list and taxonomy
fishbase.species <- rfishbase::species(server="fishbase")
fishbase.taxonomy <- rfishbase::load_taxa(server="fishbase")
fishbase.synonyms <- rfishbase::synonyms(server="fishbase")

# removed subspecies and get uniques
uniqs <- uk.list %>% 
    filter(source!="common") %>% 
    mutate(speciesName=paste(str_split_fixed(species," ",3)[,1],str_split_fixed(species," ",3)[,2])) %>% 
    distinct(speciesName)

# annotate with SpecCode
fishbase.synonyms.acc <- fishbase.synonyms %>% mutate(TaxonLevel=str_replace_all(TaxonLevel,"^species","Species")) %>% filter(Status=="accepted name" & TaxonLevel=="Species")
fishbase.synonyms.syn <- fishbase.synonyms %>% mutate(Status=str_replace_all(Status,"Synonym","synonym")) %>% filter(Status=="synonym")
# fishbase.synonyms.syn %>% pull(Status) %>% unique

# annotate with synonym status
uniqs.ann <- uniqs %>% 
    mutate(status=pull(fishbase.synonyms.acc,Status)[match(speciesName,pull(fishbase.synonyms.acc,synonym))]) %>% 
    mutate(status=if_else(is.na(status),pull(fishbase.synonyms.syn,Status)[match(speciesName,pull(fishbase.synonyms.syn,synonym))],status))

# print those without syns
drops <- uniqs.ann %>% filter(is.na(status)) %>% pull(speciesName) %>% paste(collapse=", ")
writeLines(paste0("Synonyms/accepted names for the following species could not be found in FishBase, and have been dropped: ", "\n", drops))

# drop those without status
uniqs.ann %<>% filter(!is.na(status))

# annotate with SpecCode
uniqs.ann %<>% mutate(fbSpecCode=pull(fishbase.synonyms.acc,SpecCode)[match(speciesName,pull(fishbase.synonyms.acc,synonym))]) %>% 
    mutate(fbSpecCode=if_else(is.na(fbSpecCode),pull(fishbase.synonyms.syn,SpecCode)[match(speciesName,pull(fishbase.synonyms.syn,synonym))],fbSpecCode))

# get the accepted name 
uniqs.ann %<>% mutate(validName=pull(fishbase.synonyms.acc,synonym)[match(fbSpecCode,pull(fishbase.synonyms.acc,SpecCode))])

# get synonyms
uniqs.codes <- uniqs.ann %>% pull(fbSpecCode) %>% unique()

# grab from FB - remove ones we already have, remove missapplied names and ambig synonyms, clean up any non-alphabetical chars, add valid name
fishbase.synonyms.matched <- fishbase.synonyms %>% 
    mutate(Status=str_replace_all(Status,"Synonym","synonym")) %>% 
    filter(Status=="synonym") %>%
    filter(SpecCode %in% uniqs.codes) %>% 
    filter(!synonym %in% pull(uniqs.ann,speciesName)) %>% 
    filter(!str_detect(synonym,"[^a-zA-Z\\d\\s:]")) %>% 
    mutate(validName=pull(fishbase.synonyms.acc,synonym)[match(SpecCode,pull(fishbase.synonyms.acc,SpecCode))]) %>%
    select(synonym,Status,SpecCode,validName) %>% 
    rename(speciesName=synonym,status=Status,fbSpecCode=SpecCode) 

# join with our data
spps.with.synonyms <- bind_rows(uniqs.ann,fishbase.synonyms.matched)

# keep unique rows
spps.with.synonyms %<>% distinct()

# add taxonomy and common names
spps.with.synonyms %<>% mutate(class=pull(fishbase.taxonomy,Class)[match(fbSpecCode,pull(fishbase.taxonomy,SpecCode))],
    class=pull(fishbase.taxonomy,Class)[match(fbSpecCode,pull(fishbase.taxonomy,SpecCode))],
    order=pull(fishbase.taxonomy,Order)[match(fbSpecCode,pull(fishbase.taxonomy,SpecCode))],
    family=pull(fishbase.taxonomy,Family)[match(fbSpecCode,pull(fishbase.taxonomy,SpecCode))],
    genus=pull(fishbase.taxonomy,Genus)[match(fbSpecCode,pull(fishbase.taxonomy,SpecCode))],
    commonName=pull(fishbase.species,FBname)[match(fbSpecCode,pull(fishbase.species,SpecCode))])

# make common spp table
common.spp <- uk.list %>% 
    filter(source=="common") %>% 
    mutate(fbSpecCode=pull(fishbase.synonyms.acc,SpecCode)[match(species,pull(fishbase.synonyms.acc,synonym))]) %>% 
    mutate(fbSpecCode=if_else(is.na(fbSpecCode),pull(fishbase.synonyms.syn,SpecCode)[match(species,pull(fishbase.synonyms.syn,synonym))],fbSpecCode)) %>%
    mutate(validName=pull(fishbase.species,Species)[match(fbSpecCode,pull(fishbase.species,SpecCode))])

# arrange and clean
spps.with.synonyms %<>% mutate(commonSpecies=if_else(fbSpecCode %in% pull(common.spp,fbSpecCode),TRUE,FALSE)) %>% 
    select(speciesName,status,fbSpecCode,validName,class,order,family,genus,commonName,commonSpecies) %>% 
    arrange(class,order,family,genus,validName)

# write out
write_csv(spps.with.synonyms, path="../species/uk-species-table.csv")
