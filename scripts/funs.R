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
    spp.amp <- unique(refs$sciNameValid[match(mfe$X1, refs$dbid)])
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
            paste(access.nb[a:b], collapse = ","), "&rettype=fasta&retmode=text","&", api.key)
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
                paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", 
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
