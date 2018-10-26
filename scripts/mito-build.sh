#!/usr/bin/env/sh

# check for crappy data if required (e.g. mis-specified 12S regions)
#vsearch --threads 8 --cluster_fast 12s.fas --id 0.70 --clusters ../temp/clusters/clusters.fas

GENE="12s"
GENE="coi"
FRAG="miya.primers"
FRAG="miya.noprimers"
FRAG="lerayxt.primers"
FRAG="lerayxt.noprimers"
FRAG="seamid.primers"
FRAG="seamid.noprimers"
FRAG="seashort.primers"
FRAG="seashort.noprimers"


# mafft auto align
mafft --auto --thread 8 ../hmms/mitogenome.$GENE.unaligned.fas > ../hmms/mitogenome.$GENE.aligned.fas

# now open in Geneious and cut out the Miya frag Geneious by hand
# save as e.g '../hmms/12s.miya.primers.fas'

# build the HMMs
hmmbuild ../hmms/$GENE.$FRAG.hmm ../hmms/$GENE.$FRAG.fas

# then run 'speciesDownload.R'
