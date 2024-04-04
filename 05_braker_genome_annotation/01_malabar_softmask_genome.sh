#!/bin/bash

ml bedtools/2.30.0

cd /malabar/09_Final_resources/

repeatsP0=/malabar/03_repeat_annotation/merged_repeats_unique.gff3
genomeP0=Emal_V1_24Chr_P0_rn.fasta

bedtools maskfasta -soft -fi ${genomeP0} -bed ${repeatsP0} -fo ${genomeP0/_rn/_rn_sm}
