#!/bin/bash

### Environment ###
ml bedtools/2.30.0

### Directories and Files ###
DIRin=/malabar/02_HiC_scaffolded_assembly
FILErepeats=/malabar/03_repeat_annotation/merged_repeats_unique.gff3
FILEgenome=Emal_V1_24Chr_P0_rn.fasta

### Code ###
cd ${DIRin}
bedtools maskfasta -soft -fi ${FILEgenome} -bed ${FILErepeats} -fo ${FILEgenome/.fasta/.sm.fasta}
