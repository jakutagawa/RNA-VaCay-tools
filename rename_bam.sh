#!/usr/bin/env bash
# simple script to rename all PDC objects to bams using a manifest lookup and
# index files using samtools
bam_folder="/scratch/jakutagawa/icgc/bams/normal/*"
manifest_file="/private/groups/brookslab/jakutagawa/variant_calling/download_manifests/repository_1547937442.tsv"


shopt -s nullglob
for file in $bam_folder
do
    filename=$(basename "$file")
    ext="${filename##*.}"
    # extract and check extension
    if [ -f "$file" ] && [ "$ext" != "bam" ];
    then
        object_filename=$(basename $file)
        bam_filename=$(grep $object_filename $manifest_file | cut -d$'\t' -f4)
        complete_filename=${bam_folder%?}
        complete_filename+=$bam_filename
        echo $file
        echo $complete_filename
        mv "$file" "$complete_filename"
        samtools index $complete_filename
    fi
done
shopt -u nullglob
