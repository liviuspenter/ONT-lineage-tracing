#!/bin/bash
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 2:20:00
#SBATCH -p priority
#SBATCH -o alignment_job_%A.out

# ref /n/data2/dfci/medonc/cwu/meb521/5pGT/gencode.v41.pc_canon_edit.fa.gz

echo 'mode =' $1 # map-ont sr
echo 'ref =' $2 
echo 'dir=' $3
echo 'infile=' $4
echo 'outfile=' $5

minimap2 -aY --eqx -x $1 -t 16 --secondary=no --sam-hit-only $2 $3/$4 > $3/$5.sam
samtools sort -@16 -o $3/$5.bam $3/$5.sam
samtools index -@16 $3/$5.bam
samtools view -@16 -c $3/$5.bam
