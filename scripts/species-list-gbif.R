#!/usr/bin/env Rscript

# load functions and libs
source("funs.R")

# script makes a list of all UK marine and freshwater species of fishes (inc sharks and rays).
# uses (1) gbif records; (2) fishbase; and (3) WFD transitional species list.
# only step 1 (gbif) is carried out here, fishbase and wfd list need to be done by hand.


## Get list from GBIF
# don't need to run once already done 
# just load up "../species/uk-species-list.csv" (see "scripts/species-list-clean.R")

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
gbif.sp <- gbif.db %>% mutate(species=paste(genus, specificEpithet), source=rep("gbif")) %>% select(species,source) %>% distinct(species, .keep_all=TRUE) %>% arrange(species)
#write_csv(gbif.sp, path="../temp/gbif_tmp.csv")


## Once GBIF is done, add to list from FishBase and the WFD 
# these both need to be done BY HAND from website and pdf
