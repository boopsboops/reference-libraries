#!/usr/bin/env Rscript
library("stringr")
library("seqinr")
library("parallel")
library("ape")
library("rfishbase")
library("traits")
library("rentrez")
library("bold")
library("rgbif")
library("phangorn")
library("ips")
library("stringdist")
library("tidyverse")
library("magrittr")
library("spider")
library("lubridate")
library("vroom")

# print the R session info
#sink("../temp/RsessionInfo.txt")
#print(sessionInfo())
#sink()


# function for making fasta files from tables
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")


# collapses haplotypes (from a dataframe format to a dataframe format)
# need to specify columns that contain sequence lengths, and nucleotides
# hap_collapse_df(df=mydataframe,lengthcol="lengthFrag",nuccol="nucleotidesFrag")
# add a number of each haplotype
hap_collapse_df <- function(df,lengthcol,nuccol,cores){
    odf <- df[order(df[[lengthcol]],decreasing=TRUE),]
    reps <- mcmapply(FUN=function(x) which(str_detect(string=odf[[nuccol]], pattern=x) == TRUE)[1], odf[[nuccol]], SIMPLIFY=TRUE, USE.NAMES=FALSE, mc.cores=cores)
    ind <- unique(reps)
    dat <- odf[ind,]
    dat[["nHaps"]] <- as.numeric(table(reps))
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
# results processing function to:
# remove long and short hits (20% of amplicon length)
# add taxonomy
# get best PPC per hit
# add marker and common spp etc
process_MFE <- function(mfe,primers,prefixes,references,common){
    len <- primers %>% filter(name==prefixes) %>% slice(1) %>% pull(lengthAmplicon)
    lenMax <- len+(len*0.2)
    lenMin <- len-(len*0.2)
    # filter results
    mfe %<>% filter(PPC>=0) %>% 
        filter(Size>lenMin & Size<lenMax) %>% 
        mutate(dbid=as.character(HitID)) %>%
        group_by(dbid) %>% 
        summarise(bestPPC=max(PPC)) %>% 
        ungroup() %>%
        mutate(maxPPC=max(bestPPC)) %>%
        mutate(bestPPC=(bestPPC/maxPPC)*100)
    # join
    mfe.annotated <- left_join(references,mfe)
    # choose best PPC
    mfe.annotated %<>% select(dbid,sciNameValid,class,bestPPC) %>% 
        mutate(bestPPC=if_else(is.na(bestPPC),0,bestPPC)) %>% 
        group_by(sciNameValid) %>% 
        arrange(desc(bestPPC),.by_group=TRUE) %>% 
        slice(1) %>% 
        ungroup() %>% 
        mutate(marker=prefixes) %>% 
        mutate(common=if_else(sciNameValid %in% common$sciName, "common", "rare"))
    return(mfe.annotated)
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


# function to subset the reference dataframe by marker and filter by sequence length
subset_by_marker <- function(prefix,df,thresh){
    df1 <- df %>% filter(!is.na(!!as.name(paste0("nucleotidesFrag.",prefix))))
    df1 %<>% filter(!!as.name(paste0("lengthFrag.",prefix)) >= (median(!!as.name(paste0("lengthFrag.",prefix)))*thresh))
    return(df1)
}


# removes Ns from a DNAbin list
rm_ns <-function(bin){
    bin.char <- as.character(bin)
    bin.spread <- sapply(bin.char, paste, collapse="")
    bin.rep <- sapply(bin.spread, str_replace_all, "[^actg]", "")
    bin.split <- strsplit(bin.rep, "")
    bin.bin <- as.DNAbin(bin.split)
    return(bin.bin)
}

# fun to revcomp a sequence
flip <- function(x){
    revcomp <- c2s(rev(comp(s2c(x),ambiguous=TRUE)))
    return(revcomp)}

# fun to subset a reference lib for each marker
subset_nucs <- function(pref,df){
    df %<>% rename(nucleotidesFrag=!!as.name(paste0("nucleotidesFrag.",pref)), lengthFrag=!!as.name(paste0("lengthFrag.",pref)))
    df %<>% filter(!is.na(nucleotidesFrag))
    return(df)
}

# fun to annotate a reference library table with number haplotypes per species
haps2fas <- function(df){
    df <- bind_rows(mcmapply(FUN=function(x) hap_collapse_df(df=x,lengthcol="lengthFrag",nuccol="nucleotidesFrag",cores=1), split(df,pull(df,sciNameValid)), SIMPLIFY=FALSE,mc.cores=1))
    sames <- mclapply(FUN=function(x) get_sames(df=df,ids="dbid",nucs="nucleotidesFrag",sppVec="sciNameValid",query=x), pull(df,nucleotidesFrag), mc.cores=1)
    df %<>% mutate(nMatches=sapply(sames, function(x) length(unique(x))), matchTax=sapply(sames, function(x) paste(unique(x),collapse=" | ")))
    df %<>% mutate(noms=paste(dbid,str_replace_all(sciNameValid," ","_"),nHaps,sep="|")) %>% arrange(class,order,family,genus,sciNameValid,lengthFrag,dbid)
    return(df)
}

# fun to align seqs and make a phylogentic tree
phylogenize <- function(fas,prefix,binLoc){
    fas <- ips::mafft(fas,exec="mafft",method="retree 2",maxiterate=2)
    tr <- ips::raxml(fas, file=paste0("fromR-",prefix), m="GTRCAT", f="d", p=42, exec=binLoc, N=1)
    tr <- tr$bestTree
    tmp.path <- paste0("../temp/qc_",paste(month(ymd(Sys.Date()),label=TRUE),year(ymd(Sys.Date())),sep="-"))
    dir.create(path=tmp.path)
    flist <- list.files(pattern=prefix)
    file.copy(flist, paste0(tmp.path,"/",flist))
    file.remove(flist)
    write.tree(tr,file=paste0(tmp.path,"/",prefix,".nwk"))
    return(tr)
}

# fun to plot and annotate phylogenetic trees
plot_trees <- function(tr,df,prefix){
    tr <- ape::ladderize(phangorn::midpoint(tr))
    sppv <- pull(df,sciNameValid)[match(str_split_fixed(tr$tip.label,"\\|",3)[,1],pull(df,dbid))]
    monov <- spider::monophyly(tr,sppVector=sppv)
    allmono <- monov[match(sppv, unique(sppv))]
    cols <- rep("gray20",length(tr$tip.label))
    cols[which(allmono==FALSE)] <- "hotpink"
    cols[match(df$noms[which(df$nMatches>1)], tr$tip.label)] <- "green3"
    tmp.path <- paste0("../temp/qc_",paste(month(ymd(Sys.Date()),label=TRUE),year(ymd(Sys.Date())),sep="-"))
    dfs <- df %>% summarise(nSeqs=sum(nHaps),nHaps=length(nHaps),nSpp=length(unique(sciNameValid)))
    tit <- paste0(str_replace_all(prefix,"\\.noprimers",""),"\n(n=",pull(dfs,nSeqs),", n haps=",pull(dfs,nHaps),", n spp.=",pull(dfs,nSpp),")\npink = non-monophyletic species\ngreen = shared haplotypes\nscroll down for tree ...")
    pdf(file=paste0(tmp.path,"/RAxML_bestTree.",prefix,".pdf"), width=15, height=length(tr$tip.label)/10)
    plot.phylo(tr, tip.col=cols, cex=0.5, font=1, label.offset=0.01, no.margin=TRUE)
    title(tit, line=-10)
    dev.off()
}
