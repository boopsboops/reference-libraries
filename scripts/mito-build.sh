#!/usr/bin/env/sh

GENE="12s"
GENE="16s"
GENE="coi"
GENE="cytb"
GENE="nd2"
FRAG="miya.primers"
FRAG="miya.noprimers"
FRAG="lerayxt.primers"
FRAG="lerayxt.noprimers"
FRAG="seamid.primers"
FRAG="seamid.noprimers"
FRAG="seashort.primers"
FRAG="seashort.noprimers"
FRAG="valentini.primers"
FRAG="valentini.noprimers"
FRAG="riaz.primers"
FRAG="riaz.noprimers"
FRAG="taberlet.primers"
FRAG="taberlet.noprimers"
FRAG="ward.primers"
FRAG="ward.noprimers"
FRAG="berry.noprimers"
FRAG="minamoto.noprimers"
FRAG="frag.noprimers"

# mafft auto align
mafft --auto --thread 8 ../hmms/mitogenome.$GENE.unaligned.fas > ../hmms/mitogenome.$GENE.aligned.fas

# now open in Geneious and cut out the fragment by hand
# save as e.g '../hmms/12s.miya.primers.fas'

# build the HMMs
hmmbuild ../hmms/$GENE.$FRAG.hmm ../hmms/$GENE.$FRAG.fas
