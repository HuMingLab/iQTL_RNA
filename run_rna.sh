#!/bin/bash

name=$1
CC1=$2
CC2=$3
mouse=$4
fastq1=$5
fastq2=$6
OUT=
#generate genome index

if [ ! -d "/home/leeh7/WangLab_ASE_pipeline/60_mice/test/index/${CC1}x${CC2}_index/" ]; then
STAR --runMode genomeGenerate --genomeDir /index/${CC1}x${CC2}_index/ --genomeFastaFiles /${CC1}x${CC2}_genome.fa \
 --sjdbGTFfile /${CC1}x${CC2}_anno.gtf --runThreadN 20 --sjdbOverhang 100 --genomeSAsparseD 2 --genomeSuffixLengthMax 300
else
 echo "Index already exists"
fi

#unique mapping

STAR --runThreadN 8 \
        --genomeDir  /index/${CC1}x${CC2}_index \
        --genomeLoad NoSharedMemory \
        --readFilesIn $fastq1 $fastq2 \
        --readFilesCommand zcat \
        --quantMode GeneCounts \
        --sjdbScore 2 \
        --twopassMode Basic \
        --twopass1readsN -1 \
        --outFilterMultimapNmax 1 \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        --outReadsUnmapped Fastx \
        --outFileNamePrefix ${OUT}/ase_exon_counts/${name}_


#bedtools intersect

bedtools intersect -c -a /${CC1}x${CC2}_unique_exon_annotation_ASE.BED -b ${OUT}/ase_exon_counts/${name}_Aligned.sortedByCoord.out.bam  > ${OUT}/ase_exon_counts/${name}_Counts.BED

##multimapped 

if [ ! -d "/index/${CC1}_index/" ]; then
STAR --runMode genomeGenerate --genomeDir /index/${CC1}_index/ --genomeFastaFiles /${CC1}_genome.fa \
 --sjdbGTFfile /home/leeh7/iQTL/annotation_gft/renamed/${CC1}_anno.gtf --runThreadN 20 --sjdbOverhang 100 --genomeSAsparseD 2 --genomeSuffixLengthMax 300
else
 echo "Index already exists"
fi

if [ ! -d "/index/${CC2}_index/" ]; then
STAR --runMode genomeGenerate --genomeDir /index/${CC2}_index/ --genomeFastaFiles /${CC2}_genome.fa \
 --sjdbGTFfile /renamed/${CC2}_anno.gtf --runThreadN 20 --sjdbOverhang 100 --genomeSAsparseD 2 --genomeSuffixLengthMax 300
else
 echo "Index already exists"
fi

STAR --runThreadN 12 \
        --genomeDir /index/${CC1}_index/ \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${OUT}/ase_exon_counts/${name}_Unmapped.out.mate1 ${OUT}/ase_exon_counts/${name}_Unmapped.out.mate2 \
        --quantMode GeneCounts \
        --sjdbScore 2 \
        --outFilterMultimapNmax 10 \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM Unsorted \
        --limitBAMsortRAM 5000000000 \
        --outSAMunmapped None \
        --outReadsUnmapped None \
        --outFileNamePrefix ${OUT}/unmapped_F1_total_reads/STAR_${CC1}/${CC1}_${mouse}_

STAR --runThreadN 12 \
        --genomeDir /index/${CC2}_index/ \
        --genomeLoad NoSharedMemory \
        --readFilesIn ${OUT}/ase_exon_counts/${name}_Unmapped.out.mate1 ${OUT}/ase_exon_counts/${name}_Unmapped.out.mate2 \
        --quantMode GeneCounts \
        --sjdbScore 2 \
        --outFilterMultimapNmax 10 \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM Unsorted \
        --limitBAMsortRAM 5000000000 \
        --outSAMunmapped None \
        --outReadsUnmapped None \
        --outFileNamePrefix ${OUT}/unmapped_F1_total_reads/STAR_${CC2}/${CC2}_${mouse}_

#bedtools count

samtools sort ${OUT}/unmapped_F1_total_reads/STAR_${CC1}/${CC1}_${mouse}_Aligned.out.bam -o ${OUT}/unmapped_F1_total_reads/STAR_${CC1}/sorted_${CC1}_${mouse}_Aligned.out.bam
bedtools intersect -c -a /${CC1}_unique_exon_annotation_ASE.BED -b ${OUT}/unmapped_F1_total_reads/STAR_${CC1}/sorted_${CC1}_${mouse}_Aligned.out.bam  > ${OUT}/unmapped_F1_total_reads/${CC1}_${mouse}_Counts.BED

samtools sort ${OUT}/unmapped_F1_total_reads/STAR_${CC2}/${CC2}_${mouse}_Aligned.out.bam -o ${OUT}/unmapped_F1_total_reads/STAR_${CC2}/sorted_${CC2}_${mouse}_Aligned.out.bam
bedtools intersect -c -a /${CC2}_unique_exon_annotation_ASE.BED -b ${OUT}/unmapped_F1_total_reads/STAR_${CC2}/sorted_${CC2}_${mouse}_Aligned.out.bam  > ${OUT}/unmapped_F1_total_reads/${CC2}_${mouse}_Counts.BED


python merge_read_counts.py ${OUT}/ase_exon_counts/${name}_Counts.BED ${OUT}/unmapped_F1_total_reads/${CC1}_${mouse}_Counts.BED ${OUT}/unmapped_F1_total_reads/${CC2}_${mouse}_Counts.BED ${CC1} ${CC2} ${CC1}x${CC2}_${mouse} ${OUT}/unmapped_F1_total_reads/

python merge_outputs.py ${CC1} ${CC2} ${OUT}/unmapped_F1_total_reads/${CC1}x${CC2}_${mouse}_total_counts.txt ${OUT}/unmapped_F1_total_reads/${CC1}x${CC2}_${mouse}_allelic_total_reads.txt
