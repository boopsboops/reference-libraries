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

# get uniques
uniqs <- uk.list %>% filter(source!="common") %>% select(species) %>% distinct(species)

# annotate with SpecCode
fishbase.synonyms.acc <- fishbase.synonyms %>% filter(Status=="accepted name")
fishbase.synonyms.syn <- fishbase.synonyms %>% filter(Status=="synonym" | Status=="Synonym")

# annotate with synonym status
uniqs.ann <- as_tibble(uniqs) %>% mutate(status=pull(fishbase.synonyms.acc,Status)[match(species,pull(fishbase.synonyms.acc,synonym))]) %>% 
    mutate(status=if_else(is.na(status),pull(fishbase.synonyms.syn,Status)[match(species,pull(fishbase.synonyms.syn,synonym))],status)) %>% 
    mutate(status=str_replace_all(status,"Synonym","synonym"))

# print those without syns
print("Synonyms/accepted names for the following species could not be found in FishBase. They have been dropped.")
uniqs.ann %>% filter(is.na(status)) %>% pull(species)

# drop those without status
uniqs.ann %<>% filter(!is.na(status))

# annotate with SpecCode
uniqs.ann %<>% mutate(fbSpecCode=pull(fishbase.synonyms.acc,SpecCode)[match(species,pull(fishbase.synonyms.acc,synonym))]) %>% 
    mutate(fbSpecCode=if_else(is.na(fbSpecCode),pull(fishbase.synonyms.syn,SpecCode)[match(species,pull(fishbase.synonyms.syn,synonym))],fbSpecCode))

# get the accepted name
uniqs.ann %>% mutate(validName=pull(fishbase.species,Species)[match(fbSpecCode,pull(fishbase.species,SpecCode))])
    #mutate(validName=apply(str_split_fixed(validName, " ", 3)[,1:2], 1, paste, collapse=" ")) %>% 
#filter(species!=validName) %>% print(n=Inf)

fishbase.synonyms 


# get valid names using fishbase and remove mispellings and synonyms
# takes a few minutes
val.uniqs <- validate_names(species_list=uniqs$species,server="fishbase")
val.uniqs <- unique(val.uniqs)

# check
setdiff(uniqs$species, val.uniqs)# in orig list, not in validated
setdiff(val.uniqs, uniqs$species)# in validated, not in orig
rejected.names <- c(setdiff(uniqs$species, val.uniqs),setdiff(val.uniqs, uniqs$species))

# load up the fishbase data
data(fishbase)

# make a sciname to match ours 
gen.sp <- paste(fishbase$Genus, fishbase$Species)

# match positions
fb.red <- fishbase[match(val.uniqs, gen.sp),]

# choose column want to keep and make new sciname, also sort, and tidy up 
# also edits the subspecies out 
tidy.uk.list <- fb.red %>% #
    select(Class, Order, Family, Genus, Species, SpecCode, FBname) %>% #
    mutate(sciName=paste(Genus, Species), commonSpecies=rep(FALSE), synonym=rep(FALSE), sciName=paste(str_split_fixed(sciName, " ", 3)[,1],str_split_fixed(sciName, " ", 3)[,2])) %>% #
    arrange(Class, Order, Family, Genus, Species) %>% #
    rename(class=Class, order=Order, family=Family, genus=Genus, species=Species, commonName=FBname, specCode=SpecCode) %>%
    distinct(sciName, .keep_all=TRUE)

# add the common species field
common.sp <- uk.list %>% filter(source=="common")
common.sp %>% print(n=Inf)
tidy.uk.list$commonSpecies[match(common.sp$species, tidy.uk.list$sciName)] <- TRUE

# check
glimpse(tidy.uk.list)
tidy.uk.list %>% print(n=Inf)

# function to retrieve synonyms for all species
grab_synonyms <- function(synonyms,code){
syns <- synonyms %>% filter(SpecCode==code)
return(syns)
}

# get all synonyms from rfishbase
synonyms.all <- synonyms(server="fishbase")

# run function for all spp
fish.syns <-  bind_rows(lapply(tidy.uk.list$specCode, function(x) grab_synonyms(synonyms=synonyms.all,code=x)))

# subset the synonyms and the subspecific names
fish.syns <- bind_rows((fish.syns %>% filter(Status!="accepted name")), (fish.syns %>% filter(Status=="accepted name" & TaxonLevel=="Nominotypical")))

# tidy and filter
syns.table <- fish.syns %>% 
    rename(specCode=SpecCode,sciName=synonym) %>% 
    mutate(synonym=rep(TRUE)) %>% 
    distinct(sciName, .keep_all=TRUE) %>% 
    filter(!sciName %in% tidy.uk.list$sciName) %>%
    select(specCode,synonym,sciName)

# merge with all uk species 
tidy.uk.list.syns <- bind_rows(tidy.uk.list,syns.table)

# add the rejected names from the earlier validation
tidy.uk.list.syns <- bind_rows(tidy.uk.list.syns,tibble(sciName=setdiff(rejected.names,tidy.uk.list.syns$sciName),synonym=rep(TRUE)))

# write out
write_csv(tidy.uk.list.syns, path="../species/uk-species-table.csv")
