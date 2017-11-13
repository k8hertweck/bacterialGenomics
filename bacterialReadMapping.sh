#!/bin/bash

## variant calling of bacterial genomes using read mapping

## usage: bash bacterialReadMapping.sh
# 	run from project directory
#	this script will take several minutes to run and lots of text will be printed to the screen
## dependencies (pre-installed in BioLinux)
# 	bwa: read mapping, http://bio-bwa.sourceforge.net
# 	samtools: ngs manipulation, http://www.htslib.org
#	bcfutils and vcfutils.pl: variant call analysis, https://samtools.github.io/bcftools/bcftools.html

## set variable for location of project directory (where script is executed)
PROJECT=`pwd`

## set variable for location of vcfutils.pl
VCF=/usr/share/samtools

## create directory for intermediate mapping files
echo "CREATING DIRECTORY FOR MAPPING FILES"
mkdir mapping

## index reference genome
# bwa and samtools require the index to be referenced to perform comparisons with reads
cd reference
echo "INDEXING REFERENCE"
bwa index CaceATCC824.fas
samtools faidx CaceATCC824.fas
cd $PROJECT

## read mapping
# map reads to reference: find locations in known genome where reads match
echo "MAPPING READS"
bwa mem -t 2 reference/CaceATCC824.fas data/paired*R1_001.fastq.gz data/paired*R2_001.fastq.gz > mapping/mappedReads.sam
# convert sam to sorted bam format
echo "CONVERTING SAM TO BAM"
samtools view -bS mapping/mappedReads.sam | samtools sort - mapping/mappedReads.sorted
# print simple summary statistics for read mapping
echo "SUMMARIZE READ MAPPING"
samtools flagstat mapping/mappedReads.sorted.bam > results/mappedReads.summary.txt
# add depth of coverage to summary file
echo "CALCULATING DEPTH OF COVERAGE"
samtools depth mapping/mappedReads.sorted.bam | awk '{sum+=$3} END { print "Average coverage= ",sum/NR}' >> results/reference.summary.txt

## variant calling
# find SNPs in reads relative to reference
echo "CALLING VARIANTS"
samtools mpileup -ugf reference/CaceATCC824.fas mapping/mappedReads.sorted.bam > mapping/var.raw.bcf
# filter SNPs and keep only those with substantial data
echo "FILTERING SNPS"
bcftools view mapping/var.raw.bcf | $VCF/vcfutils.pl varFilter -D100 > results/var.flt.vcf
# summarize SNPs
echo "SUMMARIZE SNPs"
echo -e "QUAL\t#non-indel\t#SNPs\t#transitions\t#joint\tts/tv\t#joint/#ref #joint/#non-indel" > results/snps.out.txt
$VCF/vcfutils.pl qstats results/var.flt.vcf >> results/snps.out.txt
