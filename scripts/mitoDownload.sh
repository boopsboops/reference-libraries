#!/usr/bin/env/sh

# download the curated mitogenome data from Miya's website
wget -v http://mitofish.aori.u-tokyo.ac.jp/files/mitogenomes.zip -P ../mitogenomes/
wget -v http://mitofish.aori.u-tokyo.ac.jp/files/mitoannotations.zip -P ../mitogenomes/

# unzip
unzip ../mitogenomes/mitogenomes.zip -d ../mitogenomes/mitogenomes/
unzip ../mitogenomes/mitoannotations.zip -d ../mitogenomes/mitoannotations/

# remove zips
rm ../mitogenomes/*.zip

# remove the 'genes' files
rm ../mitogenomes/mitogenomes/*_genes.fa

# concatentate the fasta
cat ../mitogenomes/mitogenomes/*.fa > ../mitogenomes/mitogenomes/mitogenomes.fsa

# make a file of 12S gene positions
echo -e "file\tgene\tlocation\tproduct\tgrep" > ../mitogenomes/mitoannotations/indices12S.tsv
grep "12S rRNA" ../mitogenomes/mitoannotations/*.txt >> ../mitogenomes/mitoannotations/indices12S.tsv

# open Rscript mitoExtract.R
