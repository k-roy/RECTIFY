#!/bin/bash

## required programs: trim_galore, bbmap, STAR, bwa, samfixcigar, samtools
echo "STARTING QUANTSEQ ANALYSIS PIPELINE..."

GENOME_FASTA_FILE="/Users/kevinroy/Google_Drive/yeast_genome_references/S288C_reference_sequence_R64-2-1_20150113_reformatted_chromosome_names.fasta"
GFF="/Users/kevinroy/Google_Drive/yeast_genome_references/saccharomyces_cerevisiae_annotations_GCR1_corrected.gff"

STAR_100BP_GENOME="/Users/kevinroy/Google_Drive/yeast_genome_references/STAR_indexed_R64-2-1_annotated_genome_100_bp_SJDB"
mkdir $STAR_100BP_GENOME
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $STAR_100BP_GENOME --genomeFastaFiles $GENOME_FASTA_FILE --sjdbGTFfile $GFF --sjdbOverhang 49 --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene ID

SAMPLE_NAMES="WT_37 WT_42 WT_42_IVP"

DIR="/Volumes/LaCie/pAmapping/"
FASTQ_DIR=$DIR"fastq/"
TRIMMED_DIR=$DIR"trimmed/"
mkdir $TRIMMED_DIR
ALIGNED_DIR=$DIR"aligned/"
mkdir $ALIGNED_DIR
COLLAPSED_DIR=$DIR"collapsed/"
mkdir $COLLAPSED_DIR

for sample in $SAMPLE_NAMES;
do echo "starting mapping pipeline for:" 
echo $sample

## TrimGalore on 3´ends of reads
cd $TRIMMED_DIR
trim_galore --length 15 $FASTQ_DIR$sample".fastq";
R1=$TRIMMED_DIR$sample"_trimmed.fq"

## BBMap
bbmap.sh in1=$R1 ref=$GENOME_FASTA_FILE out=$ALIGNED_DIR$sample"_bbmap.bam" local=t nhtag=t mdtag=t nhtag=t xmtag=t  amtag=t  nmtag=t  xstag=fs stoptag=t lengthtag=t idtag=t inserttag=t  scoretag=t  timetag=t boundstag=t
samtools sort $ALIGNED_DIR$sample"_bbmap.bam" -o $ALIGNED_DIR$sample"_bbmap_sorted.bam"
samtools index $ALIGNED_DIR$sample"_bbmap_sorted.bam"

## STAR
STAR --runThreadN 8 --genomeDir $STAR_100BP_GENOME --readFilesIn $R1 --outFileNamePrefix $ALIGNED_DIR$sample"_STAR_" --outSAMattributes NH nM NM MD
samtools view -bS $ALIGNED_DIR$sample"_STAR_Aligned.out.sam"  > $ALIGNED_DIR$sample"_STAR.bam"
samtools sort $ALIGNED_DIR$sample"_STAR.bam" -o $ALIGNED_DIR$sample"_STAR_sorted.bam";
java -jar /Applications/jvarkit/dist/samfixcigar.jar -r $GENOME_FASTA_FILE $ALIGNED_DIR$sample"_STAR_sorted.bam" -o $ALIGNED_DIR$sample"_STAR_sorted_reformatted_cigar.bam"
samtools index $ALIGNED_DIR$sample"_STAR_sorted_reformatted_cigar.bam";

## BWA
bwa mem $GENOME_FASTA_FILE $R1 > $ALIGNED_DIR$sample"_bwa.sam"
samtools view -bS $ALIGNED_DIR$sample"_bwa.sam"  > $ALIGNED_DIR$sample"_bwa.bam"
samtools sort $ALIGNED_DIR$sample\_bwa.bam -o $ALIGNED_DIR$sample\_bwa_sorted.bam;
java -jar /Applications/jvarkit/dist/samfixcigar.jar  -r $GENOME_FASTA_FILE $ALIGNED_DIR$sample"_bwa_sorted.bam" -o $ALIGNED_DIR$sample"_bwa_sorted_reformatted_cigar.bam"
samtools index $ALIGNED_DIR$sample"_bwa_sorted_reformatted_cigar.bam";

