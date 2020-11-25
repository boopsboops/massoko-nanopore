#!/usr/bin/env sh

# install guppy CPU
# https://community.nanoporetech.com/downloads
# follow instructions on https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revo_14dec2018
# sudo apt install ont-guppy-cpu

# set working dirs
WDIR="../temp-local-only/data/Massoko"

# 2019 runs (low quality)
RUN="3mA"
RUN="3mB"
RUN="22mA"
RUN="22mB"
RUN="22mBv2"
# 2020 runs (low quality)
RUN="3m-10"
RUN="7m-8"
RUN="12m-8"
RUN="18m-9"
RUN="22m-25"

# set quality
QUAL="lowquality"
QUAL="highquality"

# make a fastq dir
mkdir -p "$WDIR"/"$RUN"/fastq/"$QUAL"
#rm -r "$WDIR"/fastq

# check version and configs
# Nov 2020 = v4.2.2
guppy_basecaller -h
guppy_basecaller --version
guppy_basecaller --print_workflows
guppy_basecaller --print_workflows | grep "FLO-FLG001" | grep "SQK-LSK109"
# info
# https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revt_14dec2018/guppy-software-overview

# run guppy (high quality)
guppy_basecaller --input_path "$WDIR"/"$RUN"/fast5 --save_path "$WDIR"/"$RUN"/fastq/"$QUAL" --flowcell FLO-FLG001 --kit SQK-LSK109 --qscore_filtering --min_qscore 7 --records_per_fastq 4000 --compress_fastq --progress_stats_frequency 60 --num_callers 1 --cpu_threads_per_caller 7

# run guppy (low quality) see https://community.nanoporetech.com/posts/help-with-guppy-cpu-versi
guppy_basecaller --input_path "$WDIR"/"$RUN"/fast5 --save_path "$WDIR"/"$RUN"/fastq/"$QUAL" --config dna_r9.4.1_450bps_fast.cfg --qscore_filtering --min_qscore 7 --records_per_fastq 4000 --compress_fastq --progress_stats_frequency 60 --num_callers 1 --cpu_threads_per_caller 7

# info 
# R9.4.1 nanopores - 1D experiments
# are 'FLO-FLGOP1' the same as 'FLO-FLG001'? - received FLO-FLGOP1 first, then FLO-FLG001

# cat a fastq
cat "$WDIR"/"$RUN"/fastq/"$QUAL"/pass/*.fastq.gz > "$WDIR"/"$RUN"/fastq/"$QUAL"/"$RUN".fastq.gz

# write up some stats
guppy_basecaller --version > "$WDIR"/"$RUN"/fastq/"$QUAL"/stats.txt

# get n reads per sample
printf "n reads:\n" >> "$WDIR"/"$RUN"/fastq/"$QUAL"/stats.txt
wc -l "$WDIR"/"$RUN"/fastq/"$QUAL"/sequencing_summary.txt >> "$WDIR"/"$RUN"/fastq/"$QUAL"/stats.txt

# get n filtered reads
printf "\nn filtered reads:\n" >> "$WDIR"/"$RUN"/fastq/"$QUAL"/stats.txt
seqkit stats -b "$WDIR"/"$RUN"/fastq/"$QUAL"/"$RUN".fastq.gz >> "$WDIR"/"$RUN"/fastq/"$QUAL"/stats.txt

# generate md5sums of fastq files
printf "\nmd5sum:\n" >> "$WDIR"/"$RUN"/fastq/"$QUAL"/stats.txt
md5sum "$WDIR"/"$RUN"/fastq/"$QUAL"/"$RUN".fastq.gz >> "$WDIR"/"$RUN"/fastq/"$QUAL"/stats.txt

# to qc the minion run
# https://github.com/roblanf/minion_qc
Rscript ~/Software/minion_qc/MinIONQC.R -i "$WDIR"/"$RUN"/fastq/"$QUAL"/sequencing_summary.txt

# tidy up
mkdir -p "$WDIR"/"$RUN"/fastq/"$QUAL"/logs
mkdir -p "$WDIR"/"$RUN"/fastq/"$QUAL"/minion_qc
mv "$WDIR"/"$RUN"/fastq/"$QUAL"/*.log "$WDIR"/"$RUN"/fastq/"$QUAL"/logs
mv "$WDIR"/"$RUN"/fastq/"$QUAL"/*.png "$WDIR"/"$RUN"/fastq/"$QUAL"/minion_qc
