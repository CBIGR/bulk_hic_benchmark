#!/bin/bash

# Prepare epigenomic tracks (.bw format)
# For each epigenomic track: 
# --> download bam file 
# --> merge bam replicates 
# --> sort and index bam file 
# --> normalize and create .bw file
# 
# Usage example: 
# screen
# bash /epiphany/data_preparation/bam2bw.sh \
# https://www.encodeproject.org/files/ENCFF353YPB/@@download/ENCFF353YPB.bam \
# https://www.encodeproject.org/files/ENCFF677MAG/@@download/ENCFF677MAG.bam \
# GM12878_H3K36me3

BAM1=$1
BAM2=$2
DOWNLOAD_FOLDER=$3
FILE_STORE_FOLDER=$4
FILE_NAME=$5

cd "${DOWNLOAD_FOLDER}"
echo "${FILE_NAME}"
if [ -z $BAM2 ]; then
    wget "${BAM1}" -O "${FILE_STORE_FOLDER}"/"${FILE_NAME}"_merged.bam
else
    wget "${BAM1}" -O "${DOWNLOAD_FOLDER}"/"${FILE_NAME}"_rep1.bam
    wget "${BAM2}" -O "${DOWNLOAD_FOLDER}"/"${FILE_NAME}"_rep2.bam
    samtools merge -f ${FILE_STORE_FOLDER}/${FILE_NAME}_merged.bam ${DOWNLOAD_FOLDER}/${FILE_NAME}_rep1.bam ${DOWNLOAD_FOLDER}/${FILE_NAME}_rep2.bam
fi
samtools sort ${FILE_STORE_FOLDER}/${FILE_NAME}_merged.bam ${FILE_STORE_FOLDER}/${FILE_NAME}_merged_sorted
samtools index ${FILE_STORE_FOLDER}/${FILE_NAME}_merged_sorted.bam

bamCoverage --bam "${FILE_STORE_FOLDER}"/"${FILE_NAME}"_merged_sorted.bam \
    -o "${FILE_STORE_FOLDER}"/"${FILE_NAME}"_merged_sorted_RPKM.bw \
    --normalizeUsing RPKM \
    --binSize 10 \
    --ignoreForNormalization chrX

rm "${DOWNLOAD_FOLDER}"/*.bam
rm "${FILE_STORE_FOLDER}"/*.bam
rm "${FILE_STORE_FOLDER}"/*.bam.bai

