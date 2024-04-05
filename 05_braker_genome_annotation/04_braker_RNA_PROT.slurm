#!/bin/bash
#SBATCH --partition=largemem
#SBATCH -t 20-00
#SBATCH --mem=200G
#SBATCH --job-name=braker2
#SBATCH --output=04_braker_RNA_PROT.out
#SBATCH --error=04_braker_RNA_PROT.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Note ###
#Takes about 12-14 days to run

### Environment ###
source ~/.bashrc
conda activate braker2
ml genemark/4.62_lic
ml prothint/2.6.0

### Directories and Files ###
FILEgenome=Emal_V1_24Chr_P0_rn_sm.fasta #soft masked genome
DIRiso=/malabar/04_isoseq_assembly/04_read_mapping/
FILEisobam=Emal_isoforms.fastq.sorted.bam
FILErnabam=Emal_larvae_mapped_sorted.bam
FILEprotein=/databases/aa_seq_annotation/uniprot_plus_fish_proteomes.fasta

DIRout=/malabar/05_braker_annotation

### Code ###
cd ${DIRout}

#braker version 2.1.6
braker.pl	--cores=8 \
			--gff3 \
			--softmasking \
			--species=e_malabaricus_braker_RNA_PROT_all \
			--genome=${FILEgenome} \
			--bam=${DIRiso}/${FILEisobam},${FILErnabam} \
			--workingdir=emal_braker_RNA_PROT_all \
			--etpmode \
			--prot_seq=${FILEprotein}
