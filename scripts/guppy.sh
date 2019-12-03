#!/usr/bin/env sh

# install guppy CPU
# follow instructions on https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revo_14dec2018
# sudo apt install ont-guppy-cpu

# set working dirs
WDIR="../temp-local-only/data/Massoko_Shallow1/3MA/20191129_1341_MN32297_ABH252_2cedd62b"
mkdir "$WDIR"/fastq
#rm -r "$WDIR"/fastq

# check version (3.4.1+ad4f8b9) and configs
guppy_basecaller -h
guppy_basecaller --version
guppy_basecaller --print_workflows
guppy_basecaller --print_workflows | grep "FLO-FLG001" | grep "SQK-LSK109"
# info
# https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revo_14dec2018/guppy-basecaller-and-guppy-basecaller-server

# run guppy (high quality)
guppy_basecaller --input_path "$WDIR"/fast5 --save_path "$WDIR"/fastq --flowcell FLO-FLG001 --kit SQK-LSK109 --qscore_filtering --min_qscore 7 --records_per_fastq 4000 --compress_fastq --progress_stats_frequency 60 --num_callers 1 --cpu_threads_per_caller 7

# run guppy (low quality) see https://community.nanoporetech.com/posts/help-with-guppy-cpu-versi
guppy_basecaller --input_path "$WDIR"/fast5 --save_path "$WDIR"/fastq --config dna_r9.4.1_450bps_fast.cfg --qscore_filtering --min_qscore 7 --records_per_fastq 4000 --compress_fastq --progress_stats_frequency 60 --num_callers 1 --cpu_threads_per_caller 7

# info 
# R9.4.1 nanopores - 1D experiments
# are 'FLO-FLGOP1' the same as 'FLO-FLG001'? - received FLO-FLGOP1 first, then FLO-FLG001