#!/usr/bin/env Rscript
library("stringr")
library("parallel")


# collapses haplotypes (from a dataframe format to a dataframe format)
# need to specify columns that contain sequence lengths, and nucleotides
# hap_collapse_df(df=mydataframe,lengthcol="lengthFrag",nuccol="nucleotidesFrag")
hap_collapse_df <- function(df,lengthcol,nuccol){
    odf <- df[order(df[[lengthcol]],decreasing=TRUE),]
    ind <- unique(mcmapply(FUN=function(x) which(str_detect(string=odf[[nuccol]], pattern=x) == TRUE)[1], odf[[nuccol]], SIMPLIFY=TRUE, USE.NAMES=FALSE, mc.cores=8))
    dat <- odf[ind,]
    return(dat)
}


# function to get retrieve species names of sequences with an identical haplotype as your query 
# works on a dataframe
# get_sames(df=mydataframe,ids="dbid",nucs="nucleotidesFrag",sppVec="sciNameValid",query=mydataframe$nucleotidesFrag[[1]])
get_sames <- function(df,ids,nucs,sppVec,query){
    per.ind <- df[[sppVec]][str_detect(df[[nucs]], query)]
    return(per.ind)
}


# function to calculate species that drop out of a dataset after length trimming
# species_lost(df=reflib,thresh=0.5)
# threshold is a proportion of the mean sequence length
species_lost <- function(df,thresh){
    removed <- df %>% filter(lengthFrag < (mean(lengthFrag)*thresh)) %>% select(sciNameValid)
    kept <- df %>% filter(lengthFrag >= (mean(lengthFrag)*thresh)) %>% select(sciNameValid)
    tot <- setdiff(removed$sciNameValid, kept$sciNameValid)
    return(tot)
}


# function to calculate sequences removed from a dataset after length trimming
# sequences_removed(df=reflib,thresh=0.5)
# threshold is a proportion of the mean sequence length
sequences_removed <- function(df,thresh){
    removed <- df %>% filter(lengthFrag < (mean(lengthFrag)*thresh)) %>% select(dbid)
    n.removed <- length(removed$dbid)
    return(n.removed)
}
