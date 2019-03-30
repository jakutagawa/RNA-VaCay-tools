#!/usr/bin/env bash
# simple script to convert sra files to bams

sra_folder="/scratch/jakutagawa/icgc/bams/gtex/sra"
sratoolkit_folder='/private/groups/brookslab/jakutagawa/tools/sratoolkit.2.9.2-centos_linux64/bin'
star_folder='/private/groups/brookslab/bin/STAR/bin/Linux_x86_64'
star_index='/scratch/jakutagawa/reference_genomes/star_genome'
fastq_folder="/scratch/jakutagawa/icgc/bams/gtex/fastq"
bam_folder="/scratch/jakutagawa/icgc/bams/gtex/star_aligned_bams"
picard_folder="/private/groups/brookslab/bin/picard-tools-1.140"


shopt -s nullglob
for file in "$sra_folder"/*
#for dir in /scratch/jakutagawa/RNA-seq/hg38/normal/*/
do

    sra_filename=$(basename "$file")
    ext="${sra_filename##*.}"
    #sra_filename2="${file%.*}"
    sra_id="${sra_filename%.*}"
    if [ -f $file ] && [ "$ext" == "sra" ];
    then
        #echo 'first thing'
        #echo $file
        #echo $sra_filename
        #echo $sra_filename2
        fastq_file1="$fastq_folder"/"$sra_id"_1.fastq
        fastq_file2="$fastq_folder"/"$sra_id"_2.fastq
        #echo $fastq_file1
        #echo $fastq_file2
        aligned_folder="$bam_folder"/"$sra_id"/
        aligned_file="$aligned_folder"Aligned.sortedByCoord.out.bam
        log_final="$aligned_folder"Log.final.out
        echo $aligned_folder
        echo $aligned_file
        echo $sra_filename

        if [ ! -d $aligned_folder ];
        then
            mkdir $aligned_folder
        fi

        cd $aligned_folder
        bam_file="$sra_id".STAR.bam
        bam_index="$bam_file".bai
        echo $bam_file
        echo $bam_index


        #if [ ! -f $aligned_file ] | [ ! -f $bam_file ];
        if [ ! -f $log_final ];
        then
            if [ ! -f $aligned_file ];
            then

                if [ ! -f $fastq_file1 ];
                then
                    cd $fastq_folder
                    echo 'converting sra to fastq for '$sra_id
                    $sratoolkit_folder/fastq-dump --split-files $file
                    cd $aligned_folder
                else
                    echo 'fastq files already exist for '$sra_id
                fi


                echo 'aligning '$sra_id' with STAR'
                $star_folder/STAR --genomeDir $star_index \
                --readFilesIn $fastq_file1 $fastq_file2 \
                --runThreadN 8 --outFilterMultimapScoreRange 1 \
                --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 \
                --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 \
                --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 \
                --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 \
                --sjdbOverhang 100 --outSAMstrandField intronMotif\
                --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within \
                --outSAMtype BAM SortedByCoordinate --twopassMode Basic \
                --twopass1readsN -1
                echo "renaming STAR output to "$bam_file
                mv Aligned.sortedByCoord.out.bam $bam_file
            fi
        else
            if [ -f $aligned_file ];
            then
                echo "already aligned, will rename"
                mv Aligned.sortedByCoord.out.bam $bam_file
            fi
        fi

        if [ ! -f $bam_index ];
        then
            echo 'indexing '$sra_id
            samtools index $bam_file
        fi

        if [ -f $bam_file ] && [ -f $bam_index ];
        then
            echo $sra_id' successfully converted to bam'
            if [ -f $fastq_file1 ];
            then
                echo 'deleting fastq files for '$sra_id
                rm $fastq_file1
                rm $fastq_file2
            fi
        else
            echo $sra_id' not successfully converted to bam'
        fi

        picard_file="${bam_file%.*}".picard.bam


        if [ -f $bam_file ];
        then
            java -jar $picard_folder/picard.jar AddOrReplaceReadGroups \
            I=$bam_file \
            O=$picard_file \
            RGID=1 \
            RGLB=library \
            RGPL=illumina \
            RGPU=barcode \
            RGSM=$sra_id

            samtools index $picard_file




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
shopt -u nullglob
