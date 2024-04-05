#!/bin/bash
#SBATCH --partition=compute
#SBATCH -t 2:00:00
#SBATCH -c 128
#SBATCH --mem=100G
#SBATCH --job-name=STAR_indexing
#SBATCH --output=01_indexing_gff.out
#SBATCH --error=01_indexing_gff.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
ml star/2.7.9a

### Directories ###
DIRin=/malabar/02_HiC_scaffolded_assembly
FILEgenome=Emal_V1_24Chr_P0_rn.fasta
DIRgff=/03_Annotation/03_annotation_genes/03_Breaker_annotation/emal_braker_RNA_PROT_all
FILEgff=final_annot.gff3
DIRout=/malabar/07_larval_analysis/02_Read_mapping/

### Code ###
mkdir -p ${DIRout}

STAR	--runThreadN 128 \
		--runMode genomeGenerate \
		--genomeDir ${DIRout} \
		--genomeFastaFiles ${DIRgenome}/${FILEgenome} \
		--genomeSAindexNbases 13 \
		--sjdbGTFfile ${DIRgff}/${FILEgff} \
		--sjdbOverhang 149