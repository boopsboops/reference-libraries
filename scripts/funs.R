#!/usr/bin/env Rscript
library("stringr")
library("parallel")
library("ape")
library("tidyverse")
library("magrittr")
library("rfishbase")
library("traits")
library("rentrez")
library("bold")
library("rgbif")
library("phangorn")
library("ips")
library("stringdist")

## NOTE ##
# until changes come through to CRAN, need to run dev version of traits and install with sudo R in the terminal 
# either 'install_github("ropensci/traits")' or for root lib locations:
# sudo R
# library(devtools); withr::with_libpaths(new = "/usr/local/lib/R/site-library/", install_github("ropensci/traits"))

# print the R session info
#sink("../temp/RsessionInfo.txt")
#print(sessionInfo())
#sink()


# function for making fasta files from tables
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")


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
    removed <- df %>% filter(lengthFrag < (median(lengthFrag)*thresh)) %>% select(sciNameValid)
    kept <- df %>% filter(lengthFrag >= (median(lengthFrag)*thresh)) %>% select(sciNameValid)
    tot <- setdiff(removed$sciNameValid, kept$sciNameValid)
    return(tot)
}


# function to calculate sequences removed from a dataset after length trimming
# sequences_removed(df=reflib,thresh=0.5)
# threshold is a proportion of the mean sequence length
sequences_removed <- function(df,thresh){
    removed <- df %>% filter(lengthFrag < (median(lengthFrag)*thresh)) %>% select(dbid)
    n.removed <- length(removed$dbid)
    return(n.removed)
}


# function to calculate primer ID rates from MFE primer results 
primerID <- function(refs,mfe){
    spp.amp <- unique(refs$sciNameValid[which(refs$dbid %in% mfe$X1)])
    spp.tot <- unique(refs$sciNameValid)
    num <- length(spp.amp)/length(spp.tot)
    names(num) <- paste0(prefix, ", ", length(spp.tot), " total species")
    return(num)
}


# `read.GenBank` function modified from ape package, but to now include an api key 
read_GenBank <- function (access.nb, seq.names = access.nb, species.names = TRUE, 
    gene.names = FALSE, as.character = FALSE, api.key) 
{
    N <- length(access.nb)
    a <- 1L
    b <- if (N > 400) 
        400L
    else N
    fl <- tempfile()
    repeat {
        URL <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
            paste(access.nb[a:b], collapse = ","), "&rettype=fasta&retmode=text&api_key=", api.key)
        X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
        cat(X, sep = "\n", file = fl, append = TRUE)
        if (b == N) 
            break
        a <- b + 1L
        b <- b + 400L
        if (b > N) 
            b <- N
    }
    res <- read.FASTA(fl)
    if (is.null(res)) 
        return(NULL)
    attr(res, "description") <- names(res)
    if (length(access.nb) != length(res)) {
        names(res) <- gsub("\\..*$", "", names(res))
        failed <- paste(access.nb[!access.nb %in% names(res)], 
            collapse = ", ")
        warning(paste0("cannot get the following sequences:\n", 
            failed))
    }
    else names(res) <- access.nb
    if (as.character) 
        res <- as.character(res)
    if (species.names) {
        a <- 1L
        b <- if (N > 400) 
            400L
        else N
        sp <- character(0)
        repeat {
            URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
                paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text&api_key=", api.key,
                sep = "")
            X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE, 
                n = -1)
            sp <- c(sp, gsub(" +ORGANISM +", "", grep("ORGANISM", 
                X, value = TRUE)))
            if (b == N) 
                break
            a <- b + 1L
            b <- b + 400L
            if (b > N) 
                b <- N
        }
        attr(res, "species") <- gsub(" ", "_", sp)
    }
    if (gene.names) 
        warning("you used 'gene.names = TRUE': this option is obsolete; please update your code.")
    res
}


# R script to run a hidden markov model on a sequence
# need to have "hmmer" and "biosquid" installed 
# if not, run 'sudo apt install hmmer biosquid'
# also needs the ape package
# requires a tempfile directory (e.g. "../temp")
# requires an infile in fasta format to be in the same dir as the tempfiles (e.g. "myfile.fas")
# requires the name of the hmm you want to use (e.g. "12s.miya.noprimers.hmm")
# requires a prefix for the hmmer output (e.g. "analysis1")
# assumes the hidden markov model is located in ../hmms 
# returns a DNAbin object of the sequences matched by hmmer 

run_hmmer3 <- function(dir, infile, hmm, prefix, evalue, coords){#
    string.hmmer <- paste0("nhmmer -E ", evalue, " --incE ", evalue, " --dfamtblout ", dir, "/", prefix, ".hmmer.tbl ", "../hmms/", prefix, ".hmm ", dir, "/", infile)
    system(command=string.hmmer, ignore.stdout=TRUE)
    hmm.tbl <- read_delim(file=paste0(dir, "/", prefix, ".hmmer.tbl"), delim=" ", col_names=FALSE, trim_ws=TRUE, progress=FALSE, comment="#", col_types=cols(), guess_max=100000)
    names(hmm.tbl) <- c("targetName","acc","queryName","bits","eValue","bias","hmmStart","hmmEnd","strand","aliStart","aliEnd","envStart","envEnd","sqLen","descriptionTarget")
    hmm.tbl %<>% filter(strand=="+") %>% distinct(targetName, .keep_all=TRUE) %>% mutate(coords=paste(envStart,envEnd,sep=":"))
    mtdna <- read.FASTA(file=paste0(dir,"/",infile))
    mtdna.sub <- as.character(mtdna[match(hmm.tbl$targetName,names(mtdna))])
    if(coords=="env"){
    mtdna.sub.coords <- as.DNAbin(mapply(function(x,y,z) x[y:z], mtdna.sub, hmm.tbl$envStart, hmm.tbl$envEnd, SIMPLIFY=TRUE, USE.NAMES=TRUE))
    } else if(coords=="ali"){
    mtdna.sub.coords <- as.DNAbin(mapply(function(x,y,z) x[y:z], mtdna.sub, hmm.tbl$aliStart, hmm.tbl$aliEnd, SIMPLIFY=TRUE, USE.NAMES=TRUE))
    } else {
    stop("Please provide 'env' or 'ali' as arguments to coords")
    }
    return(mtdna.sub.coords)
}
