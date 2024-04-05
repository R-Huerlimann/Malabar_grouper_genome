#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-12
#SBATCH -c 5
#SBATCH --mem=10G
#SBATCH --job-name=03_mapQC
#SBATCH --output=03_mapping_QC-%j.out
#SBATCH --error=03_mapping_QC-%j.err
#SBATCH --mail-user=roger.huerlimann@oist.jp
#SBATCH --mail-type=END,FAIL
#SBATCH --array=0-23%24

### Modules ###
module load samtools
module load rseqc
module load multiqc

### Directories ###
DIRin=/malabar/07_larval_analysis/02_Read_mapping/
array=(${DIRin}/REMD*)
DIRbed=/malabar/05_braker_annotation/emal_braker_RNA_PROT_all
FILEbed=final_annot.bed


### Code ###

echo "changing directory to ${array[${SLURM_ARRAY_TASK_ID}]}"
cd ${array[${SLURM_ARRAY_TASK_ID}]}

name=${array[${SLURM_ARRAY_TASK_ID}]};name=${name##*/}
echo "Processing ${name}"

file=${name}"Aligned.sortedByCoord.out.bam"

echo ""
echo "indexing bam file for ${name}"
samtools index ${file}
echo ""
echo "running sammtools flagstat on ${name}"
samtools flagstat ${file} > ${name}_flagstat.txt

echo ""
echo "running bam_stat on ${name}"
bam_stat.py -i ${file} > ${name}_bam_stat.txt

echo ""
echo "running read_distribution on ${name}"
read_distribution.py -r ${bedFile} -i ${file} > ${name}_read_distribution.txt

echo ""
echo "running geneBody_coverage on ${name}"
geneBody_coverage.py -i ${file} -r ${bedFile} -o ${name}_geneBodyCoverage

echo ""
echo "running tin on ${name}"
tinn.py -i ${file} -r ${bedFile}



#when finished, run the below in the 07_larval_analysis folder to get a summary of all QC results
#multiqc -ip --filename multiQC_mapping .
