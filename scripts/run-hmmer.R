#!/usr/bin/env Rscript
# R script to run a hidden markov model on a sequence

# need to have "hmmer" and "biosquid" installed 
# if not, run 'sudo apt install hmmer biosquid'
# also needs the ape package

# requires a tempfile directory (e.g. "../temp")
# requires an infile in fasta format to be in the same dir as the tempfiles (e.g. "myfile.fas")
# requires the name of the hmm you want to use (e.g. "12s.miya.noprimers.hmm")
# requires a prefix for the hmmer output (e.g. "analysis1")
# requires a 

# assumes the hidden markov model is located in ../data 
# currently it's called '../data/miya.hmm'

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

### OLD func using the aligned output
#run_hmmer <- function(dir, infile, hmm, prefix, evalue){#
#    string.hmmer <- paste0("nhmmer -E ", evalue, " --incE ", evalue, " -A", " ", dir, "/", prefix, ".hmmer.stk", " ", "--notextw", " ", "--tblout", " ", dir, "/", prefix, ".hmmer.tbl", " ", "-o", " ", dir, "/", prefix, ".hmmer.out", " ", "../hmms/", prefix, ".hmm", " ", dir, "/", infile)
#    string.format <- paste0("sreformat fasta", " ", dir, "/", prefix, ".hmmer.stk", " > ", dir, "/", prefix, ".hmmer.fas")
#    system(command=string.hmmer)
#    system(command=string.format)
#    nuc <- read.dna(file=paste0(dir, "/", prefix, ".hmmer.fas"), format="fasta")
#    names(nuc) <- str_split_fixed(names(nuc), " \\[subseq from\\] ", 2)[,2]
#    return(nuc)
#}#
