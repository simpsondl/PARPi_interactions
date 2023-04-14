#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=8:00:00
#SBATCH --mem=50G
#SBATCH --job-name="twostep"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out


#### Call script with                 
#### sbatch --export=workdir=WORKDIR,ref=REF,r1=R1,r2=R2,out=OUT bwa_protospacer.s

module purge
module load htslib/1.9
module load samtools

cd /Genomics/grid/users/ds65/adamsonlab/data/HTS_data/DS/HTS_JL019/align/

#Get sample info from filenames
int1=${in1##*/}
int2=${int1#demultiplexed_}
int3=${int2%.fq.gz}
repl=${int3%%_*}
read=${int3##*_}
sample=${out%%_*}

echo $sample

#Align reads to constant region
bowtie2 -p 8 -5 19 --ff -x /Genomics/adamsonlab/data/HTS_data/HTS_JL011/hts_jl011/bowtie_index/cr \
-1 $in1 \
-2 $in2 |\
samtools view -bS - > $out

#Extract read ids to keep
samtools view -F 0xc -f 0x40 -q 5 $out | awk '$3=="CR2" && $7=="CR3"' | cut -f1 > ${sample}_${repl}_crfilt_keep.txt

#Filter reads
~/tools/seqtk/seqtk subseq $in1 ${sample}_${repl}_crfilt_keep.txt > ${sample}_${repl}_r1_crfilt.fq
~/tools/seqtk/seqtk subseq $in2 ${sample}_${repl}_crfilt_keep.txt > ${sample}_${repl}_r2_crfilt.fq

gzip ${sample}_${repl}_r1_crfilt.fq 
gzip ${sample}_${repl}_r2_crfilt.fq

#Align filtered reads to protospacers
bowtie2 -p 8 -3 21 --ff -x /Genomics/adamsonlab/data/HTS_data/HTS_JL018/ref/aunip_lib_reference \
-1 ${sample}_${repl}_r1_crfilt.fq.gz \
-2 ${sample}_${repl}_r2_crfilt.fq.gz |\
samtools view -bS - > ${sample}_${repl}_ps.bam

samtools view -F 0xc -f 0x40 -q 5 ${sample}_${repl}_ps.bam |\
cut -f3,7 | sort | uniq -c > ${sample}_${repl}_counts.txt
