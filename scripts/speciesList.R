#!/usr/bin/env Rscript

# load libs
require("tidyverse")
require("rgbif")
require("rfishbase")
require("taxize")
require("stringr")
require("parallel")

# script makes a list of all UK marine and freshwater species of fishes (inc sharks and rays)
# uses (1) gbif records; (2) fishbase; and (3) WFD transitional species list.


## Get list from GBIF
# Don't need to run once already done (done 20/12/17)

# get taxa names for sharks and fish
keyfish <- name_backbone(name='Actinopterygii')$classKey
keyshark <- name_backbone(name='Elasmobranchii')$classKey

# get max counts for records (only preserved specimen records with locality data, i.e. with vouchers) 
maxfish <- occ_count(taxonKey=keyfish, country="GB", georeferenced=TRUE, basisOfRecord="PRESERVED_SPECIMEN")
maxshark <- occ_count(taxonKey=keyshark, country="GB", georeferenced=TRUE, basisOfRecord="PRESERVED_SPECIMEN")

# download the records (only preserved specimen records with locality data, i.e. with vouchers) 
fishdb <- occ_search(taxonKey=keyfish, country="GB", hasCoordinate=TRUE, basisOfRecord="PRESERVED_SPECIMEN", return="data", fields="all", limit=maxfish)
sharkdb <- occ_search(taxonKey=keyshark, country="GB", hasCoordinate=TRUE, basisOfRecord="PRESERVED_SPECIMEN", return="data", fields="all", limit=maxshark)

# dataset dimensions
dim(sharkdb)
dim(fishdb)

# merge the dbs
gbif.db <- dplyr::full_join(fishdb, sharkdb)

# removed those with no species level ID and those with non-euro locations (overseas territories and mistakes) and extant from palao collections
gbif.db <- gbif.db %>% filter(taxonRank == "SPECIES") %>% filter(!(decimalLatitude < 40 | decimalLongitude < -30)) %>% filter(collectionCode!="PAL")

# get the names and write out to temp 
gbif.sp <- gbif.db %>% mutate(species=paste(genus, specificEpithet)) %>% select(species) %>% distinct(species, .keep_all=TRUE) %>% arrange(species)
#write_csv(gbif.sp, path="../temp/temp/gbif_tmp.csv")


## Once GBIF is done add to list from FishBase and the WFD 
# these both need to be done BY HAND from website and pdf

# load up the combined species list
uk.list <- read_csv(file="../data/uk_species_list.csv")

# get uniques
uniqs <- uk.list %>% filter(source!="common") %>% select(species) %>% distinct(species)

# get valid names using fishbase and remove mispellings and synonyms
# takes a few minutes
val.uniqs <- validate_names(species_list=uniqs$species)
val.uniqs <- unique(val.uniqs)
# check for probs
warnings()

# check
setdiff(uniqs$species, val.uniqs)
setdiff(val.uniqs, uniqs$species)

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

# now get synonyms for all species
fish.syns <- mclapply(tidy.uk.list$specCode, FUN=function(x) rfishbase::synonyms(x,limit=100), mc.cores=8)

# join and clean
syns.table <- bind_rows(fish.syns) %>% 
    mutate(genSpp=paste(SynGenus, SynSpecies), sciName=paste(str_split_fixed(genSpp, " ", 3)[,1], str_split_fixed(genSpp, " ", 3)[,2]), synonym=rep(TRUE)) %>% 
    distinct(sciName, .keep_all=TRUE) %>% 
    filter(!sciName %in% tidy.uk.list$sciName)

# reduce
syns.table <- syns.table %>% select(SpecCode,SynGenus,SynSpecies,sciName,synonym) %>% rename(specCode=SpecCode,genus=SynGenus,species=SynSpecies)

# merge with all uk species
tidy.uk.list.syns <- bind_rows(tidy.uk.list,syns.table)

# write out
write_csv(tidy.uk.list.syns, path="../data/uk_species_table.csv")
