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

run_hmmer <- function(dir, infile, hmm, prefix, evalue){#
    string.hmmer <- paste0("nhmmer -E ", evalue, " -A", " ", dir, "/", prefix, ".hmmer.stk", " ", "--notextw", " ", "--tblout", " ", dir, "/", prefix, ".hmmer.tbl", " ", "-o", " ", dir, "/", prefix, ".hmmer.out", " ", "../hmms/", prefix, ".hmm", " ", dir, "/", infile)
    string.format <- paste0("sreformat fasta", " ", dir, "/", prefix, ".hmmer.stk", " > ", dir, "/", prefix, ".hmmer.fas")
    system(command=string.hmmer)
    system(command=string.format)
    nuc <- read.dna(file=paste0(dir, "/", prefix, ".hmmer.fas"), format="fasta")
    names(nuc) <- sapply(strsplit(names(nuc), split="/"), function(x) x[1])
    return(nuc)
}#
