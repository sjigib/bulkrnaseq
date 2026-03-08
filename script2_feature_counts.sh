#!/bin/bash

#SBATCH --job-name=feature
#SBATCH --output=feature.out
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --cpus-per-task=30
#SBATCH --nodelist=node12
#SBATCH --qos=common
#SBATCH --account=common

SECONDS=0

for bam_file in /lustre/user/rna-seq_pipeline/star/*Aligned.sortedByCoord.out.bam; do 
    # Get the sample name from the BAM file name
     sample_name=$(basename "$bam_file" | sed 's/Aligned.sortedByCoord.out.bam//')

    # Run featureCounts command
    featureCounts -p -T 8 -t exon -g gene_id -a /lustre/user/human_ref_genome/gencode.v44.basic.annotation.gtf -o /lustre/user/rna-seq_pipeline/star/feature_counts/"$sample_name"_counts.txt "$bam_file"
done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "featurecount generated successfully!"












