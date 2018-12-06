#!/usr/bin/env Rscript

# this script cleans and annotates the uk fish species list using fishbase

# load functions and libs
source("funs.R")

# load up the combined species list (generated 20/12/2017)
uk.list <- read_csv(file="../species/uk-species-list.csv")

# get uniques
uniqs <- uk.list %>% filter(source!="common") %>% select(species) %>% distinct(species)

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
