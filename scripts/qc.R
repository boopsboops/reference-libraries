#!/usr/bin/env Rscript
# script to quality control the reference libraries and identify erroneous sequences.

# load functions and libs
source("funs.R")

# load up the reference lib
reflib <- read_csv(file="../temp/libraries_29-10-18/uk-fishes.12s.miya.noprimers.csv",guess_max=100000)

# plot the length distributions
reflib %>% ggplot(aes(lengthFrag)) + geom_histogram(binwidth=1)

# make a dataframe of trimming stats
trimming.df <- data_frame(percMed=seq(0,1,by=0.05), 
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

# the suspect sequences are added to the exclusions file
exclusions <- read_csv(file="../references/exclusions.csv")

# format some temp names
# make col with seq duplicates
sames <- mclapply(FUN=function(x) get_sames(df=reflib,ids="dbid",nucs="nucleotidesFrag",sppVec="sciNameValid",query=x), reflib$nucleotidesFrag, mc.cores=8)
reflib <- reflib %>% mutate(nMatches=sapply(sames, function(x) length(unique(x))), matchTax=sapply(sames, function(x) paste(unique(x),collapse=" | ")))

reflib.tmp <- reflib %>% mutate(noms=paste(dbid,str_replace_all(sciNameValid," ","_"),sep="|")) %>% arrange(class,order,family,genus,sciNameValid,lengthFrag,dbid)
# make fasta 
reflib.fas <- tab2fas(df=reflib.tmp,seqcol="nucleotidesFrag",namecol="noms")
#write.FASTA(reflib.fas, file="../temp/temp/riaz.fas")# write if needed

# align the sequences with MAFFT (need to have exe in PATH)
sam <- mafft(reflib.fas,path="mafft",method="retree 1")

# make a quick ML tree with RAxML ((need to have exe on your system))
# ignore the 'cannot open file' error
# takes about 12 hours for full COI!
raxml(sam, m="GTRCAT", f="d", p=111118, exec="~/Software/standard-RAxML/raxmlHPC-AVX", N=1)

# read in the tree
rax.tr <- read.tree("RAxML_bestTree.fromR")
rax.tr <- midpoint(ladderize(rax.tr))

# copy or delete the log info raxml created - and move the tree file elsewhere
file.rename(from="RAxML_bestTree.fromR", to="../../temp/primer-faceoff/raxml/RAxML_bestTree.miya.nwk")
file.remove(dir(path=".", pattern="fromR"))

# color tips
cols <- vector("character", length(rax.tr$tip.label))
cols[match(reflib.tmp$noms[which(reflib.tmp$nMatches>1)], rax.tr$tip.label)] <- "blue"
cols[cols!="blue"] <- "black"

# plot PDF
# adjust margins
pdf(file="../../temp/primer-faceoff/raxml/miya.all.pdf", width=15, height=190)
plot.phylo(rax.tr, tip.col=cols, cex=0.5, font=1, label.offset=0.01, no.margin=TRUE)
dev.off()


# to make a Levenstein distance tree
seqs.char <- as.character(lapply(as.character(reflib.fas), str_flatten))
names(seqs.char) <- names(lapply(as.character(reflib.fas), str_flatten))

# construct matrix - takes a LONG time on COI data
smat <- stringdistmatrix(seqs.char, method="lv", useNames="names")
rax.tr <- midpoint(ladderize(nj(smat)))
# go to prev code to format and print tree