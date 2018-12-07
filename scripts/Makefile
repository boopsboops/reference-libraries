# makefile for generating species reference library scripts from GenBank searches and tissues collected
# Rupert Collins, Aug 2018
# http://swcarpentry.github.io/make-novice/02-makefiles/

# format
# output : dependences
# commands (preceded by a tab)

# to execute, simply run "make" in the terminal from the scripts dir
 
../docs/reference-library-tables.md : reference-library-tables.Rmd ../species/tissues.csv ../references/uk-fish-references.csv.gz
	Rscript -e "require('rmarkdown'); rmarkdown::render(input='reference-library-tables.Rmd', output_dir='../docs/')"
