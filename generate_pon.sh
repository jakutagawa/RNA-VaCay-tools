#!/usr/bin/env bash
# simple script to convert sra files to bams

#sra_folder="/scratch/jakutagawa/icgc/bams/gtex/sra"
sratoolkit_folder='/private/groups/brookslab/jakutagawa/tools/sratoolkit.2.9.2-centos_linux64/bin'
#star_folder='/private/groups/brookslab/bin/STAR/bin/Linux_x86_64'
star_index='/scratch/jakutagawa/reference_genomes/star_genome'
#fastq_folder="/scratch/jakutagawa/icgc/bams/gtex/fastq"
ref_fasta="/private/groups/brookslab/reference_indices/hs37"
bam_folder="/scratch/jakutagawa/icgc/bams/gtex/star_aligned_bams"
picard_folder="/private/groups/brookslab/bin/picard-tools-1.140"
gatk_folder="/private/groups/brookslab/jakutagawa/tools/gatk-4.0.11.0/gatk"
pon_vcf_folder='/scratch/jakutagawa/icgc/bams/gtex/pon'
vcf_file_list=$pon_vcf_folder/normals_for_pon_vcf.args


shopt -s nullglob
#for file in "$sra_folder"/*
for dir in $bam_folder/*
do
    #echo $dir
    for file in $dir/*
    do
        #echo $file

        filename=$(basename "$file")

        ext="${filename##*.}"
        ext2="${filename#*.}"
        #sra_filename2="${file%.*}"

        if [ -f $file ] && [ "$ext2" == "STAR.bam" ];
        then
            sra_id="${filename%.STAR.bam}"
            echo $filename
            echo $sra_id
            echo $ext2
            #echo 'first thing'
            #echo $file
            #echo $sra_filename
            #echo $sra_filename2

            #echo $fastq_file1
            #echo $fastq_file2


            bam_file="$sra_id".STAR.bam
            bam_index="$bam_file".bai
            echo $bam_file
            echo $bam_index


            picard_file="${file%.*}".picard.bam
            picard_index=$picard_file.bai
            split_picard="${file%.*}".picard.split.bam

            vcf_out=$pon_vcf_folder/"$sra_id".vcf.gz

            echo $picard_file
            echo $picard_index
            echo $vcf_out

            if [ ! -f $vcf_out ];
            then
                if [ ! -f $picard_file ];
                then
                    echo 'replacing header on '$sra_id
                    java -jar $picard_folder/picard.jar AddOrReplaceReadGroups \
                    I=$file \
                    O=$picard_file \
                    RGLB=library \
                    RGPL=illumina \
                    RGPU=barcode \
                    RGSM=$sra_id
                else
                    echo 'picard file exists'
                fi

                if [ ! -f $picard_index ];
                then
                    echo 'indexing '$sra_id
                    samtools index $picard_file
                fi

                if [ ! -f $split_picard ];
                then
                    echo 'splitting file with SplitNCigarReads'
                    $gatk_folder SplitNCigarReads \
                    -R $ref_fasta/hs37d5.fa \
                    -I $picard_file \
                    -O $split_picard
                fi

                echo 'running '$sra_id' on Mutect2'
                #echo $picard_file
                #echo $ref_fasta/hs37d5.fa
                #echo $vcf_out
                $gatk_folder Mutect2 \
                -R $ref_fasta/hs37d5.fa \
                -I $split_picard \
                -tumor $sra_id \
                -O $vcf_out
                echo $vcf_out >> $vcf_file_list

            else
                echo $sra_id' vcf already exists'
                echo $vcf_out >> $vcf_file_list
            fi

            #if [ -f $picard_file ];
            #then
            #    rm $picard_file
            #fi
            #if [ -f $picard_index ];
            #then
            #    rm $picard_index
            #fi

        #for dir2 in "$dir"/*

        #    for file in "$dir2"/*
        #    do
        #        if [[ -f $file ]]
        #        then
        #            echo $file
        #            #do_something_with "$file"
        #        fi
        #    done
        #done
        fi
    done
done
shopt -u nullglob
$gatk_folder CreateSomaticPanelOfNormals \
    -vcfs $vcf_file_list \
    -O $pon_vcf_folder/pon.vcf.gz
#cat $vcf_file_list
