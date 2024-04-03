#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-4
#SBATCH -c 128
#SBATCH --mem=150G
#SBATCH --job-name=BUSCO
#SBATCH --output=02_QC_BUSCO.out
#SBATCH --error=02_QC_BUSC.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp


### Environment ###
module load BUSCO/4.1.2
export AUGUSTUS_CONFIG_PATH="/BUSCO_DB/config/"

### Directories ###
DIRin=/malabar/01_assembly/
FILEgenome=/malabar/01_assembly/19-final/final.p_ctg.fasta

### BUSCO QC ###
busco   -m genome \
        -i ${FILEgenome} \
        -c 128 \
        -l /busco_downloads/lineages/actinopterygii_odb10/ \
        --out_${DIRin} \
        -o BUSCO_Results \
        -f
