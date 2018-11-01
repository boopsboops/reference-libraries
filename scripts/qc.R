#!/usr/bin/env Rscript
# script to qulaity control the reference libraries and identify erroneous sequences.

# load funs and libs
library("tidyverse")
source("funs.R")

# load up the reference lib
reflib <- read_csv(file="../temp/libraries_29-10-18/uk-fishes.12s.miya.noprimers.csv",guess_max=100000)

# plot the length distributions
reflib %>% ggplot(aes(lengthFrag)) + geom_histogram(binwidth=1)

# make a dataframe of trimming stats
trimming.df <- data_frame(percMean=seq(0,1,by=0.05), 
    bp=sapply(seq(0,1,by=0.05), function(x) round(mean(reflib$lengthFrag)*x)),
    speciesLost=sapply(seq(0,1,by=0.05), function(x) length(species_lost(df=reflib,thresh=x))),
    seqsRemoved=sapply(seq(0,1,by=0.05), function(x) sequences_removed(df=reflib,thresh=x)))

# print trim results
trimming.df %>% print(n=Inf)

# look at lost species at each value
sapply(seq(0,1,by=0.05), function(x) species_lost(df=reflib,thresh=x))
