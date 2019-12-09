#!/usr/bin/env sh

# install kraken2
git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh bin

# build a standard db
kraken2-build --standard --threads 8 --db ../temp-local-only/refdb/standard

# download a custom db (bacteria)
export PATH=~/Software/kraken2/bin:$PATH
kraken2-build --download-taxonomy --skip-maps --db ../temp-local-only/refdb/bacteria
kraken2-build --download-library bacteria --db ../temp-local-only/refdb/bacteria
kraken2-build --build --threads 8 --db ../temp-local-only/refdb/bacteria
# clean
kraken2-build --clean --db ../temp-local-only/refdb/bacteria

# premade
# https://lomanlab.github.io/mockcommunity/mc_databases.html
# https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads
mkdir ../temp-local-only/refdb/kraken2-microbial-fatfree/
cd ../temp-local-only/refdb/kraken2-microbial-fatfree/
wget -c https://refdb.s3.climb.ac.uk/kraken2-microbial/hash.k2d
wget https://refdb.s3.climb.ac.uk/kraken2-microbial/opts.k2d
wget https://refdb.s3.climb.ac.uk/kraken2-microbial/taxo.k2d

# run kraken2
WDIR="../temp-local-only/data/Massoko"
RUN="3mA"
RUN="3mB"
RUN="22mA"
RUN="22mB"
RUN="22mBv2"

# run
kraken2 --threads 8 --gzip-compressed --output ../temp-local-only/results/"$RUN".out.tsv --report ../temp-local-only/results/"$RUN".report.tsv --db ../temp-local-only/refdb/kraken2-microbial-fatfree "$WDIR"/"$RUN"/fastq/"$RUN".fastq.gz

# add sample number by sed for summary report
sed "s/^/$RUN\t/g" ../temp-local-only/results/"$RUN".report.tsv | sed "s/[\t][ ]*/\t/g" > ../temp-local-only/results/"$RUN".report.clean.tsv
cat ../temp-local-only/results/*.report.clean.tsv > ../temp-local-only/results/kraken.summary.tsv

# for by read results
sed "s/^/$RUN\t/g" ../temp-local-only/results/"$RUN".out.tsv > ../temp-local-only/results/"$RUN".out.clean.tsv
cat ../temp-local-only/results/*.out.clean.tsv > ../temp-local-only/results/kraken.results-by-read.tsv


# viz
# https://genomics.sschmeier.com/ngs-taxonomic-investigation/index.html#kraken2