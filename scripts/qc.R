#!/usr/bin/env Rscript
# script to quality control the reference libraries and identify erroneous sequences.

# load functions and libs
rm(list=ls())
# load functions
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/funs.R")
# load reference lib
source("https://raw.githubusercontent.com/boopsboops/reference-libraries/master/scripts/references-load.R")
# loads objects: uk.species.table, uk.species.table.common, reflib.orig

# set cores
cores <- 8

# make a copy so don't have to keep reloading the reflib
reflib <- reflib.orig


# list primer prefixes
prefix <- "coi.lerayxt.noprimers"
prefix <- "coi.seamid.noprimers"
prefix <- "coi.seashort.noprimers"
prefix <- "coi.ward.noprimers"
prefix <- "12s.miya.noprimers"
prefix <- "12s.riaz.noprimers"
prefix <- "12s.valentini.noprimers"
prefix <- "12s.taberlet.noprimers"
prefix <- "16s.berry.noprimers"
prefix <- "cytb.minamoto.noprimers"

# subset by primer, renaming the fragment columns
reflib %<>% rename(nucleotidesFrag=!!as.name(paste0("nucleotidesFrag.",prefix)), lengthFrag=!!as.name(paste0("lengthFrag.",prefix)))

# subset DNAs
reflib %<>% filter(!is.na(nucleotidesFrag))


## Stats on seq lengths
# plot the length distributions
reflib %>% ggplot(aes(lengthFrag)) + geom_histogram(binwidth=1)

# make a dataframe of trimming stats
trimming.df <- tibble(percMed=seq(0,1,by=0.05), 
    bp=sapply(seq(0,1,by=0.05), function(x) round(median(reflib$lengthFrag)*x)),
    speciesLost=sapply(seq(0,1,by=0.05), function(x) length(species_lost(df=reflib,thresh=x))),
    seqsRemoved=sapply(seq(0,1,by=0.05), function(x) sequences_removed(df=reflib,thresh=x)))

# print trim results
trimming.df %>% print(n=Inf)

# look at lost species at each value
lsp <- sapply(seq(0,1,by=0.05), function(x) species_lost(df=reflib,thresh=x))
names(lsp) <- paste0((seq(0,1,by=0.05)*100), "% of median seq len")
print(lsp)


## make some trees to identify erroneous sequences
# criteria for pruning:

# haplotype identical to one of another species, or clustered within obvious group
# number of sequences and number of independent studies suggest mislabelling
# lack of information in BOLD or GenBank record 
# phylogenetically implausible (clear sample mixup)

# collapse haps by spp
reflib.red <- bind_rows(mcmapply(FUN=function(x) hap_collapse_df(df=x,lengthcol="lengthFrag",nuccol="nucleotidesFrag"), split(reflib,reflib$sciNameValid), SIMPLIFY=FALSE, mc.cores=cores))

# format some temp names
# make col with seq duplicates
sames <- mclapply(FUN=function(x) get_sames(df=reflib.red,ids="dbid",nucs="nucleotidesFrag",sppVec="sciNameValid",query=x), reflib.red$nucleotidesFrag, mc.cores=cores)
reflib.red %<>% mutate(nMatches=sapply(sames, function(x) length(unique(x))), matchTax=sapply(sames, function(x) paste(unique(x),collapse=" | ")))

reflib.tmp <- reflib.red %>% mutate(noms=paste(dbid,str_replace_all(sciNameValid," ","_"),nHaps,sep="|")) %>% arrange(class,order,family,genus,sciNameValid,lengthFrag,dbid)
# make fasta 
reflib.fas <- tab2fas(df=reflib.tmp,seqcol="nucleotidesFrag",namecol="noms")
#write.FASTA(reflib.fas, file="../temp/temp/berry.fas")# write if needed
#add some unassigned seqs
#reflib.fas <- c(reflib.fas,tab2fas(df=bad.matches, seqcol="dnas", namecol="qseqid"))

## <----
## read in the generated 12S data from our tissue sample collection to add to and QC also 
## skip if you just dealing with GenBank records

# filter 12S
refs.tissue <- read.FASTA(file="../../SeaDNA/data/reference-library.fasta")
refs.tissue <- refs.tissue[grep("12S",names(refs.tissue))]
names(refs.tissue) <- str_replace_all(names(refs.tissue), "12S\\|", "")
# generate names
tissues.master <- read_csv(file="../species/tissues-master.csv")
tissues.master %<>% mutate(sciName=paste(genus,specificEpithet,sep="_"))
names(refs.tissue) <- paste(names(refs.tissue), tissues.master$sciName[match(names(refs.tissue), tissues.master$otherCatalogNumbers)], sep="|")
# sort
refs.tissue.sorted <- refs.tissue[order(str_split_fixed(names(refs.tissue),"\\|",2)[,2])]

# list of species missing from reflib
got <- str_replace_all(names(refs.tissue.sorted), "\\|.*", "")

spp.got <- tissues.master %>% filter(otherCatalogNumbers %in% got) %>% pull(sciName) %>% unique() %>% str_replace_all("_"," ")

tissues.master %>% filter(phylum!="Arthropoda",!is.na(locality),!is.na(specificEpithet)) %>% 
    filter(!sciName %in% spp.got) %>% 
    select(otherCatalogNumbers,sciName,locality) #%>% #print(n=Inf)
    #write_csv(path="../../SeaDNA/temp/reference-library/outstanding_species.csv")

# write out fas quick
write.FASTA(refs.tissue.sorted, file="../temp/temp/12Sreflib.fas")

# get just the miya frag with hmm
dat.frag <- run_hmmer3(dir="../temp/temp", infile="12Sreflib.fas", prefix="12s.miya.noprimers", evalue="10", coords="env")

# join with master reflib
reflib.fas <- c(reflib.fas,dat.frag)
## < ---


# align the sequences with MAFFT (need to have exe in PATH)
sam <- ips::mafft(reflib.fas,exec="mafft",method="retree 1")
# write out matrix to check alignment
# write.FASTA(sam,file="../temp/temp/12S.matrix.fas")
# sam <- sam[,299:542]
# as.character(sam[1,])

# make a quick ML tree with RAxML ((need to have exe on your system))
# ignore the 'cannot open file' error
# takes several hours for COI full
rax.tr <- ips::raxml(sam, m="GTRCAT", f="d", p=42, exec="~/Software/standard-RAxML/raxmlHPC-AVX", N=1)

# ladderise
rax.tr <- midpoint(ladderize(rax.tr$bestTree))

# clean up after raxml
write.tree(rax.tr,file=paste0("../../SeaDNA/temp/primer-faceoff/raxml/RAxML_bestTree.",prefix,".nwk"))
file.remove(dir(path=".", pattern="fromR"))

# color tips
cols <- vector("character", length(rax.tr$tip.label))
cols[match(reflib.tmp$noms[which(reflib.tmp$nMatches>1)], rax.tr$tip.label)] <- "blue"
cols[cols!="blue"] <- "black"
#cols[grep("^otu",rax.tr$tip.label)] <- "red"

# plot PDF
# adjust margins
pdf(file=paste0("../../SeaDNA/temp/primer-faceoff/raxml/RAxML_bestTree.",prefix,".pdf"), width=15, height=70)#900 for COI
plot.phylo(rax.tr, tip.col=cols, cex=0.5, font=1, label.offset=0.01, no.margin=TRUE)
dev.off()


# to make a Levenstein distance tree
#seqs.char <- as.character(lapply(as.character(reflib.fas), str_flatten))
#names(seqs.char) <- names(lapply(as.character(reflib.fas), str_flatten))

# construct matrix - takes a LONG time on COI data
#smat <- stringdistmatrix(seqs.char, method="lv", useNames="names")
#rax.tr <- midpoint(ladderize(nj(smat)))
# go to prev code to format and print tree