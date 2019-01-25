#!/usr/bin/env/sh

# create a temp directory
mkdir ../temp

# download the curated mitogenome data from Miya's website
wget -v http://mitofish.aori.u-tokyo.ac.jp/files/mitogenomes.zip -P ../temp/
wget -v http://mitofish.aori.u-tokyo.ac.jp/files/mitoannotations.zip -P ../temp/

# unzip
unzip ../temp/mitogenomes.zip -d ../temp/mitogenomes/
unzip ../temp/mitoannotations.zip -d ../temp/mitoannotations/

# remove zips
rm ../temp/*.zip

# remove the 'genes' files
rm ../temp/mitogenomes/*_genes.fa

# concatentate the fasta
cat ../temp/mitogenomes/*.fa > ../temp/mitogenomes/mitogenomes.fsa

# make a file of 12S gene positions
echo -e "file\tgene\tlocation\tproduct\tgrep" > ../temp/mitoannotations/indices12S.tsv
grep "12S rRNA" ../temp/mitoannotations/*.txt >> ../temp/mitoannotations/indices12S.tsv

# make a file of COI gene positions
echo -e "file\tcds\tlocation\tcodon\tposition" > ../temp/mitoannotations/indicesCOI.tsv
grep -B 1 -A 1 "COI$" ../temp/mitoannotations/*.txt | grep "codon_start" >> ../temp/mitoannotations/indicesCOI.tsv

echo -e "file\tgene\tlocation\tproduct\tgrep" > ../temp/mitoannotations/indices16S.tsv
grep "16S rRNA" ../temp/mitoannotations/*.txt >> ../temp/mitoannotations/indices16S.tsv

echo -e "file\tgene\tlocation\tproduct\tgrep" > ../temp/mitoannotations/indicesCYTB.tsv
grep -B 1 -A 1 "Cyt b$" ../temp/mitoannotations/*.txt | grep "codon_start" >> ../temp/mitoannotations/indicesCYTB.tsv

# open Rscript mitoExtract.R
