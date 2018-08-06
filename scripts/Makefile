# makefile for generating species reference library scripts from GenBank searches and tissues collected
# Rupert Collins, Aug 2018
# http://swcarpentry.github.io/make-novice/02-makefiles/

# format
# output : dependences
# commands (preceded by a tab)

# to execute, simply run "make" in the terminal from the scripts dir
 
../docs/reference_library_tables.md : reference_library_tables.Rmd ../species/tissues.csv ../references/uk-fishes-miya-12s-noprimers.csv
	Rscript -e "require('rmarkdown'); rmarkdown::render(input='reference_library_tables.Rmd', output_dir='../docs/')"
