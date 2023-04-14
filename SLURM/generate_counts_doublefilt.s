#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=150G
#SBATCH --job-name="count"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out


#### Call script with                 
#### sbatch --export=workdir=WORKDIR,ref=REF,r1=R1,r2=R2,out=OUT bwa_protospacer.s

module purge
module load samtools 

cd $indir

echo Starting in $indir

samtools view -F 0xc -f 0x80 $in |\
awk '$17=="NM:i:0"||$17=="NM:i:1"' |\
cut -f1 > ${prefix}_filt_keep.txt

java -Xmx64g -jar ~/tools/picard.jar FilterSamReads \
I=${in} \
O=${prefix}_doublefilt.bam \
READ_LIST_FILE=${prefix}_filt_keep.txt \
FILTER=includeReadList

samtools view -F 0xc -f 0x40 ${prefix}_doublefilt.bam |\
awk '$17=="NM:i:0"||$17=="NM:i:1"' | cut -f3,7 |\
sort -T /scratch/tmp/ds65/tmp/ | uniq -c > ${prefix}_doublefilt_counts.txt
