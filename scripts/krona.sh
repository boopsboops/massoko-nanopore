#!/usr/bin/env sh

# install Krona and NCBI dbs
git clone https://github.com/marbl/Krona.git
cd Krona/KronaTools
sudo ./install.pl
mkdir taxonomy
./updateTaxonomy.sh
./updateAccessions.sh

# go to R script 'kraken-plot.R' and export krona data

# run krona
ktImportTaxonomy -o ../temp/krona-taxonomy.html ../temp-local-only/results/krona-22m-combined.tsv ../temp-local-only/results/krona-3m-combined.tsv