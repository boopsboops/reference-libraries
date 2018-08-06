#!/usr/bin/env Rscript

# libs
library(tidyverse)
library(ape)

# read index file and fasta file
tab.indices <- read_tsv(file="../mitogenomes/mitoannotations/indices12S.tsv")
all.mito <- read.dna("../mitogenomes/mitogenomes/mitogenomes.fsa", format="fasta", as.character=TRUE)

# strip out the file path
tab.indices <- tab.indices %>% mutate(file=str_replace_all(file,"../mitogenomes/mitoannotations/",""))

# break up the file name
split.names <- str_split_fixed(tab.indices$file,"_",3)

# clean and make new columns
tab.indices <- tab.indices %>% mutate(accession=paste0(split.names[,1],"_",split.names[,2]), 
    species=str_replace_all(split.names[,3], ".txt:", ""), 
    species=str_replace_all(species, "_", " "), 
    location=str_replace_all(location,"\\.\\.",":"))

# do the same for the fasta using only the accession number
names(all.mito) <- str_replace_all(str_split_fixed(names(all.mito),"\\|",7)[,4], "\\.[0-9]", "")

# check names are in the same order
summary(names(all.mito) == tab.indices$accession)

# add the DNA into the table
#tab.indices %>% mutate(nucleotides=match())

# index the locations to get 12S
unaligned.12s <- as.DNAbin(mapply(function(dna,locations){dna[eval(parse(text=locations))]}, dna=all.mito, locations=tab.indices$location))

# write out
write.dna(unaligned.12s, file="../hmms/mitogenome.12s.unaligned.fas", format="fasta", colw=999999)

# now align using 'mitoAlign.sh'

# primers
MiFishUF <- "GTCGGTAAAACTCGTGCCAGC"
MiFishUR <- "CAAACTGGGATTAGATACCCCACTATG"#rev="CATAGTGGGGTATCTAATCCCAGTTTG"

