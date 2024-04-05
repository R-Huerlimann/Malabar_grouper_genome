#!/bin/bash
#SBATCH -p compute
#SBATCH -t 3-0
#SBATCH -c 5
#SBATCH --mem=20G
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Note ###
# The final libraries were split cross two NovaSeq lanes.
# Here I merge the outputs from the two lanes into one

### Data and folders ###
DIRinA=/malabar/00_Raw_data/RNAseq_Larvae/ID240_A
DIRinB=/malabar/00_Raw_data/RNAseq_Larvae/ID240_B

DIRout=/malabar/00_Raw_data/RNAseq_Larvae/ID240_AB

### Code ###
mdkir ${DIRout}

for file in ${DIRinA}/*gz; do
	zcat ${file} ${file/A/B} | gzip  > ${file/A/AB}
done
