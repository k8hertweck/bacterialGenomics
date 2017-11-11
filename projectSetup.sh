#!/bin/bash

## installing software and downloading data

# usage: bash setup.sh URL_R1 URL_R2

# assign variable for location of project directory
PROJECT=`pwd`

# download and install trimmomatic
cd # change to home directory (appropriate for software installation)
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
# if you already have trimmomatic installed, select A when prompted for an answer
unzip Trimmomatic-0.36.zip
# you can test to see if trimmomatic installed correctly by running the following:
# java -jar ~/Trimmomatic-0.36/trimmomatic-0.36.jar -version
# if the software is correctly installed, "0.36" should be printed to the screen

# setup data directory
cd $PROJECT
mkdir data

# download data from URL
wget data/$1
wget data/$2

## check data download: view first few lines of each fastq file
#for x in data/*.fastq.gz
#	do
#		echo $x
#		zcat $x | head
#done
