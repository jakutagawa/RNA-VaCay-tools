#!/usr/bin/env bash
# script to rename all call objects to bams using a manifest lookup and
# index files using samtools
bam_folder="/scratch/jakutagawa/icgc/bams/normal/*"
manifest_file="/private/groups/brookslab/jakutagawa/variant_calling/download_manifests/repository_1547937442.tsv"
output_bam_folder="/scratch/jakutagawa/icgc/bams/synthetic/"
project_codes_file="/private/groups/brookslab/jakutagawa/variant_calling/reference_files/all_project_codes.txt"
PCAWG_maf_folder="/private/groups/brookslab/PCAWG/Oct2016_Freeze/*"
custom_brca_file="/private/groups/brookslab/jakutagawa/variant_calling/reference_files/histo_Breast-CombinedCa_October_2016_whitelist_2583.snv_mnv_indel.maf"
randomized_mutation_folder="/private/groups/brookslab/jakutagawa/variant_calling/synthetic_mutation_lists/"
seed=1

shopt -s nullglob
for bam_file in $bam_folder
do
    filename=$(basename "$bam_file")
    ext="${filename##*.}"
    # extract and check extension
    #echo $bam_filename
    if [ -f "$bam_file" ] && [ "$ext" == "bam" ];
    then
        object_filename=$(basename $bam_file)
        project=$(grep $object_filename $manifest_file | cut -d$'\t' -f10)
        complete_filename=${bam_folder%?}
        complete_filename+=$bam_filename
        echo $project
        project_code=$(grep $project $project_codes_file | cut -d$'\t' -f1)
        echo $project_code
        #mv "$bam_file" "$complete_filename"
        #samtools index $complete_filename
        if [ "$project" == "BRCA-US" ];
        then
            maf_filename=$custom_brca_file
        else
            maf_filename='/private/groups/brookslab/PCAWG/Oct2016_Freeze/histo_'
            maf_filename+=$project_code
            maf_filename+='_October_2016_whitelist_2583.snv_mnv_indel.maf'
        fi
        echo $maf_filename
        echo $object_filename

        randomized_mutation_list_file=$randomized_mutation_folder
        randomized_mutation_list_file+=$(basename ${object_filename%.*})
        randomized_mutation_list_file+='.random_snp_only.txt'
        #echo $randomized_mutation_list_file

        #echo $seed
        echo $bam_file
        output_bam=$output_bam_folder
        output_bam+=${object_filename%.*}
        output_bam+='.with_variants.bam'
        echo $output_bam

        if [ ! -f "$randomized_mutation_list_file" ];
        then
            echo 'would run randomization'
            python generate_random_mutations.py -t SNP -mf $maf_filename -n 300 -s $seed -vc 5'UTR 3'UTR Missense_Mutation Nonsense_Mutation Nonstop_Mutation Silent Start_Codon_SNP > $randomized_mutation_list_file
        else
            echo 'randomized file exists'
        fi

        if [ ! -f "$output_bam" ];
        then
            echo 'begin mutation spike in '$bam_file
            python samvar.py -vf $randomized_mutation_list_file -ib $bam_file -ob $output_bam
            #samtools index $output_bam
        else
            echo 'spiked bam exists'
        fi
        ((seed++))
    fi
done
shopt -u nullglob
