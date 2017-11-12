#!/bin/bash

## installing software and downloading data

# usage: bash setup.sh URL_R1 URL_R2

# assign variable for location of project directory
PROJECT=`pwd`

# download and install trimmomatic
cd # change to home directory (appropriate for software installation)
echo " SETTING UP TRIMMOMATIC"
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
# if you already have trimmomatic installed, select A when prompted for an answer
unzip Trimmomatic-0.36.zip
# you can test to see if trimmomatic installed correctly by running the following:
# java -jar ~/Trimmomatic-0.36/trimmomatic-0.36.jar -version
# if the software is correctly installed, "0.36" should be printed to the screen

# setup data directory
echo "SETTING UP DATA DIRECTORY"
cd $PROJECT
mkdir data

# download data from URL
echo "DOWNLOADING DATA"
wget $1
wget $2
echo "MOVING DATA"
mv *.fastq.gz data/

## check data download: view first few lines of each fastq file
#for x in data/*.fastq.gz
#	do
#		echo $x
#		zcat $x | head
#done
