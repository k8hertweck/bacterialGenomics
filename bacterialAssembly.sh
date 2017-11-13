#!bin/bash

## genomics for bacterial sequences

## usage: bash bacterialAssembly.sh FILENAME_R1_001.fastq.gz FILENAME_R2_001.fastq.gz
# 	run from project directory
# 	additional explanations included in comments below
## dependencies (installed with projectSetup.sh):
# 	Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic 
## dependencies (pre-installed in BioLinux)
# 	fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# 	velvet: https://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf
#	blast: 

## set variables for read data files
R1=$1
R2=$2

## set variable for location of project directory (where script is executed)
PROJECT=`pwd`

## set variable for Trimmomatic (trimmomatic installed using setup.sh)
TRIMMOMATIC=~/Trimmomatic-0.36

# set up project directory for analysis
cd $PROJECT
mkdir results results/fastqc

## move to data directory
cd $PROJECT/data

## check quality of raw sequence data
echo "CHECKING SEQUENCE QUALITY"
# run fastqc to summarize each fastq file
fastqc $R1 $R2
mv *.html *.zip $PROJECT/results
# see more options using fastqc -h
# summarize quality assessment of raw data
cd $PROJECT/results
echo -n > data_summary.txt
for x in Cace*.zip
	do
		echo $x >> data_summary.txt
		# unzip file
		unzip $x
		# print fastqc summaries to screen
		cat *fastqc/summary.txt >> data_summary.txt
		# remove old files
		rm -rf *fastqc
done
cd $PROJECT/data

## filter and trim sequences by quality
# slashes at end of lines allow entering multi-line commands more easily
# if entering below chunk manually, copy and paste each line without slash
java -jar $TRIMMOMATIC/trimmomatic-0.36.jar PE -threads 2 -phred33 $R1 $R2 paired-$R1 unpaired-$R1 paired-$R2 unpaired-$R2 ILLUMINACLIP:$TRIMMOMATIC/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:15 TRAILING:15 HEADCROP:10 MINLEN:50

## recheck quality of raw sequence data
fastqc paired-*.fastq.gz
mv *.html *.zip $PROJECT/results
# summarize quality assessment of trimmed data (paired files only)
cd $PROJECT/results
for x in paired-*.zip
	do
		echo $x >> data_summary.txt
		unzip $x
		cat *fastqc/summary.txt >> data_summary.txt
		rm -rf *fastqc
done

## move back to main project directory
cd $PROJECT

## assemble sequence data
# run velveth, specifying kmer and read information
# find more info about command with: velveth -h
velveth velvetOut 31 -shortPaired -fastq -separate data/paired-$R1 data/paired-$R2
# run velvetg, specifying location of output from velveth
# find more info about command with: velvetg -h
velvetg velvetOut
# copy resulting contigs to results
cp velvetOut/contigs.fa results/contigs.fa

## create blast database from velvet assembly
makeblastdb -in velvetOut/contigs.fa -dbtype nucl

## search (blast) for genes of interest in assembly contigs
# genes from chromosome: HBD, CRT, BCD, ETFA, ETFB
# genes from plasmid: ADHE, THIL
cd reference
for genes in chromosome_genes.fas plasmid_genes.fas
	do
		blastn -query $genes -db $PROJECT/velvetOut/contigs.fa -outfmt 7 -out $PROJECT/results/BLAST-$genes
done

