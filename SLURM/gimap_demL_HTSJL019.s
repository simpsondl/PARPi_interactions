#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=8:00:00
#SBATCH --mem=50G
#SBATCH --job-name="deml"
#SBATCH --output=/Genomics/grid/users/ds65/logs/%x-%j.out


#### Call script with                 
#### sbatch --export=workdir=WORKDIR,ref=REF,r1=R1,r2=R2,out=OUT bwa_protospacer.s

module purge

~/tools/deML/src/deML -i /Genomics/adamsonlab/data/HTS_data/DS/HTS_JL019/demultiplexed/sample_index.txt \
-f /Genomics/adamsonlab/data/HTS_data/DS/HTS_JL019/fastq/4939__HTS_JL019-for-96-cycles-000000000-DKGBM_1_Read_1_passed_filter.fastq.gz \
-r /Genomics/adamsonlab/data/HTS_data/DS/HTS_JL019/fastq/4939__HTS_JL019-for-96-cycles-000000000-DKGBM_1_Read_2_Index_Read_passed_filter.fastq.gz \
-if1 /Genomics/adamsonlab/data/HTS_data/DS/HTS_JL019/fastq/4939__HTS_JL019-for-96-cycles-000000000-DKGBM_1_Read_3_Index_Read_passed_filter.fastq.gz \
-if2 /Genomics/adamsonlab/data/HTS_data/DS/HTS_JL019/fastq/4939__HTS_JL019-for-96-cycles-000000000-DKGBM_1_Read_4_passed_filter.fastq.gz \
-o /Genomics/adamsonlab/data/HTS_data/DS/HTS_JL019/demultiplexed/demultiplexed \
-s /Genomics/adamsonlab/data/HTS_data/DS/HTS_JL019/demultiplexed/summary.txt \
-e /Genomics/adamsonlab/data/HTS_data/DS/HTS_JL019/demultiplexed/error.txt
