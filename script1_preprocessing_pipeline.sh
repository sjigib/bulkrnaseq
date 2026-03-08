#!/bin/bash

#SBATCH --job-name=pipeline
#SBATCH --output=rnaseq.out
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --cpus-per-task=40
#SBATCH --nodelist=node12
#SBATCH --qos=common
#SBATCH --account=common

SECONDS=0

ids=("SRR15852393" "SRR15852394" "SRR15852395" "SRR15852396" "SRR15852397" "SRR15852398" "SRR15852399" "SRR15852400" "SRR15852401" "SRR15852402" "SRR15852423" "SRR15852424" "SRR15852425" "SRR15852426" "SRR15852427" "SRR15852428" "SRR15852429" "SRR15852430" "SRR15852431" "SRR15852432")

for id in "${ids[@]}"; do
# Downloading files Using SRAToolkit 
fasterq-dump --split-files "$id" --outdir /lustre/user/rna-seq_pipeline/fastq_seq

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "Raw fastq file downloaded successfully!"

# Zipping the fastq files
gzip /lustre/user/rna-seq_pipeline/fastq_seq/"$id"_1.fastq
gzip /lustre/user/rna-seq_pipeline/fastq_seq/"$id"_2.fastq

duration=$SECONDS 
echo "$(($duration/60)) minutes and $(($duration%60)) seconds."
echo "Files zipped !"

# Checking fastqc for files before trimming
fastqc -t 40 /lustre/user/rna-seq_pipeline/fastq_seq/"$id"_1.fastq.gz -o /lustre/user/rna-seq_pipeline/fastqc
fastqc -t 40 /lustre/user/rna-seq_pipeline/fastq_seq/"$id"_2.fastq.gz -o /lustre/user/rna-seq_pipeline/fastqc

duration=$SECONDS 
echo "$(($duration/60)) minutes and $(($duration%60)) seconds."
echo "Quality check before trimming completed !"

# Trimming adapters using Trimmomatic
# Specify the input directory containing your paired-end read files
input_dir="/lustre/user/rna-seq_pipeline/fastq_seq"

# Specify the output directory for the processed reads
output_dir="/lustre/user/rna-seq_pipeline/trimmed"

# Run Timmomatic with the input and output file paths for the current sample
trimmomatic PE -threads 8 "$input_dir/$id"_1.fastq.gz "$input_dir/$id"_2.fastq.gz "$output_dir/$id"_1.paired.trimmed.fastq.gz "$output_dir/$id"_1.unpaired.trimmed.fastq.gz "$output_dir/$id"_2.paired.trimmed.fastq.gz "$output_dir/"$id_2.unpaired.trimmed.fastq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

duration=$SECONDS
echo "$((${duration} / 60)) minutes and $((${duration} % 60)) seconds elapsed."
echo "Trimming done successfully!"

# Checking fastqc for files after trimming
fastqc -t 40 /lustre/user/rna-seq_pipeline/trimmed/"$id"_1.paired.trimmed.fastq.gz -o /lustre/user/rna-seq_pipeline/fastqc2
fastqc -t 40 /lustre/user/rna-seq_pipeline/trimmed/"$id"_2.paired.trimmed.fastq.gz -o /lustre/user/rna-seq_pipeline/fastqc2

duration=$SECONDS 
echo "$(($duration/60)) minutes and $(($duration%60)) seconds."
echo "Quality check after trimming completed !"

# Aligning files using STAR
# Set the path to the reference genome index
REF_GENOME_INDEX="/lustre/user/human_ref_genome/ref_genome"

# Set the input directory containing your paired-end read files
INPUT_DIR="/lustre/user/rna-seq_pipeline/trimmed"

# Set the output directory for the aligned reads
OUTPUT_DIR="/lustre/user/rna-seq_pipeline/star"

# Run STAR alignment for the current sample
STAR --runThreadN 20 --genomeDir "$REF_GENOME_INDEX" --readFilesCommand gunzip -c --readFilesIn "$INPUT_DIR/$id"_1.paired.trimmed.fastq.gz "$INPUT_DIR/$id"_2.paired.trimmed.fastq.gz --outFileNamePrefix "$OUTPUT_DIR/$id"_ --outSAMtype BAM SortedByCoordinate  

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "alignment run successfully!"

done