# iQTL pipline for RNA-seq
## Allelic RNA-seq mapping for CCGG.

Find the Hi-C pipeline [here](ttps://github.com/lindsayhrlee/iQTL_HiC).
Find the ATAC-seq pipeline [here](https://github.com/lindsayhrlee/iQTL_ATAC).

### Requirements

- STAR v2.7.9a
- Bedtools v2.29.1
- Samtools v1.9
- Python v3.6.8
    * Numpy v1.19.5
    * Pandas v1.1.5
 
### Required inputs
- fasta file (with STAR index)
- gtf file
- annotation bed file
- The fastq file should be named CC1xCC2_mousename_R1_*.fastq.gz


### Running the pipeline
This pipeline can be run simply by:

```
./run_rna.sh [Mouse] [CC1] [CC2] [Name] [fastq R1] [fastq R2]
```

 The required input variables are:
1.	Mouse = The mouse ID (CC1xCC2_Name)
2.	CC1 = Maternal Genome Name
3.	CC2 = paternal Genome Name 
4.	Name = Mouse name
5.	fastq R1 = FASTQ file for R1
6.	fastq R2 = FASTQ file for R2

Here are the detailed steps to this pipeline: (Be sure to fill out all the locations for index, gff, fasta, and annotation files)

Step 1: Indexing pseudo-genome fasta file using STAR
	Input: pseudo-genome fasta file, annotation gtf file
	Output: A folder that contains STAR indexing information
```
STAR --runMode genomeGenerate --genomeDir /index/${CC1}x${CC2}_index/ --genomeFastaFiles  CC1}x${CC2}_genome.fa \
 --sjdbGTFfile ${ CC1}x${CC2}_anno.gtf --runThreadN 20 --sjdbOverhang 100 --genomeSAsparseD 2 --genomeSuffixLengthMax 300
```
Step 2: Unique mapping via STAR
	Input: Index folder from step 1, R1/R2 fastq files
	Output: a BAM file, unmapped reads files that will be used for multimapped mapping step
```
STAR --runThreadN 8 \
        --genomeDir  /home/leeh7/WangLab_ASE_pipeline/60_mice/test/index/${CC1}x${CC2}_index \
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
```
Step 3: Getting the read counts via bedtools
	Input: mapped BAM file from step 2, exon annotation file
	Output: a BED file that contains read counts
```
bedtools intersect -c -a /home/leeh7/iQTL/gff/name_added_gtf/unique/${CC1}x${CC2}_unique_exon_annotation_ASE.BED -b ${OUT}/ase_exon_counts/${name}_Aligned.sortedByCoord.out.bam  > ${OUT}/ase_exon_counts/${name}_Counts.BED
```
Step 4: indexing maternal and paternal fasta files separately via STAR
	Input: maternal/paternal fasta files, respective annotation file 
	Output: A folder that contains STAR indexing information
```
STAR --runMode genomeGenerate --genomeDir ${CC1}_index/ --genomeFastaFiles ${CC1}_genome.fa \
--sjdbGTFfile ${CC1}_anno.gtf --runThreadN 20 --sjdbOverhang 100 --genomeSAsparseD 2 --genomeSuffixLengthMax 30

STAR --runMode genomeGenerate --genomeDir ${CC2}_index/ --genomeFastaFiles ${CC2}_genome.fa \
 --sjdbGTFfile  ${CC2}_anno.gtf --runThreadN 20 --sjdbOverhang 100 --genomeSAsparseD 2 --genomeSuffixLengthMax 300
```
Step 5: Multimapping maternal and paternal genomes separately 
	Input: Maternal and Paternal Index folder, R1/R2 fastq files
	Output: a BAM file
```
STAR --runThreadN 12 \
        --genomeDir index/${ CC1}_index/ \
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
        --genomeDir  /index/${CC2}_index/ \
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
```
Step 6: sorting the BAM files and counting the reads
	Input: BAM files from step 5, exon annotation file
	Output: a BED file that contains read counts
```
samtools sort ${OUT}/unmapped_F1_total_reads/STAR_${CC1}/${CC1}_${mouse}_Aligned.out.bam -o ${OUT}/unmapped_F1_total_reads/STAR_${CC1}/sorted_${CC1}_${mouse}_Aligned.out.bam
bedtools intersect -c -a /${CC1}_unique_exon_annotation_ASE.BED -b ${OUT}/unmapped_F1_total_reads/STAR_${CC1}/sorted_${CC1}_${mouse}_Aligned.out.bam  > ${OUT}/unmapped_F1_total_reads/${CC1}_${mouse}_Counts.BED

samtools sort ${OUT}/unmapped_F1_total_reads/STAR_${CC2}/${CC2}_${mouse}_Aligned.out.bam -o ${OUT}/unmapped_F1_total_reads/STAR_${CC2}/sorted_${CC2}_${mouse}_Aligned.out.bam
bedtools intersect -c -a /${CC2}_unique_exon_annotation_ASE.BED -b ${OUT}/unmapped_F1_total_reads/STAR_${CC2}/sorted_${CC2}_${mouse}_Aligned.out.bam  > ${OUT}/unmapped_F1_total_reads/${CC2}_${mouse}_Counts.BED
```
Step 7: Merging all the read counts
	Input: Read counts bed files from unique and multimapped 
	Output: a merged total counts file
```
python merge_read_counts.py ${OUT}/ase_exon_counts/${name}_Counts.BED ${OUT}/unmapped_F1_total_reads/${CC1}_${mouse}_Counts.BED ${OUT}/unmapped_F1_total_reads/${CC2}_${mouse}_Counts.BED ${CC1} ${CC2} ${CC1}x${CC2}_${mouse} ${OUT}/unmapped_F1_total_reads/
```
Step 8: Including the gene information in the merged total counts file
	Input: the merged total counts file from step 7
   Output: final output that contains the gene information along with the read counts
```
python merge_outputs.py ${CC1} ${CC2} ${OUT}/unmapped_F1_total_reads/${CC1}x${CC2}_${mouse}_total_counts.txt ${OUT}/unmapped_F1_total_reads/${CC1}x${CC2}_${mouse}_allelic_total_reads.txt
```


### Contact Us

For any questions regarding this software, contact Ming Hu (hum@ccf.org) or Lindsay Lee (leeh7@ccf.org).
