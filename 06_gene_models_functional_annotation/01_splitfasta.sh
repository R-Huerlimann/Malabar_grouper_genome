#!/bin/bash

### Functions ###
split_fasta(){
	input=$1
	species_in=$2
	cat $input | bioawk -c 'fastx' -v species_in=$2 '{ if( (NR-1)%1000==0 ){file=sprintf("%s%d.fa",species_in,(NR-1));} printf(">%s\t%s\n%s\n",$name,$comment,$seq) >> file}'
}

### Directories and Files ###
DIRin=/malabar/05_braker_annotation/emal_braker_RNA_PROT_all
FILEaa=augustus.hints_anysupport_nosupportwithdomaindiamond_longest_isoform.aa
name=malabar
DIRout=/malabar/06_gene_models_functional_annotation

### code ###
cd ${DIRout}
split_fasta ${inDir}/${FILEaa} ${name}
