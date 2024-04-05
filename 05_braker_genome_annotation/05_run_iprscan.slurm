#!/bin/bash
#SBATCH --partition=compute
#SBATCH -t 72:00:00
#SBATCH --mem=200G
#SBATCH -c 8
#SBATCH --job-name=05_iprscan
#SBATCH --output=05_run_iprscan-%j.out
#SBATCH --error=05_run_iprscan-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roger.huerlimann@oist.jp

### Environment ###
ml interproscan/5.48.83.0

### Directories and Files ###
DIRin=/malabar/03_Annotation/04_braker_all_202109/emal_braker_RNA_PROT_all/
FILEaa=augustus.hints.aa

### Code ###
cd ${DIRin}

# Remove * at the end of sequence first
cat ${FILEaa} | sed 's|*||' > ${FILEaa/.aa/_wo_asterisk.aa/}

interproscan.sh -i ${FILEaa/.aa/_wo_asterisk.aa/} -o iprscan_result.txt -f tsv -appl Pfam -cpu 8
