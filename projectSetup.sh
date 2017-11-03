#!/bin/bash

## installing software and downloading data
# usage: bash setup.sh URL_R1 URL_R2

PROJECT=`pwd`

# download and install trimmomatic
cd
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
# if you already have trimmomatic installed, select A to overwrite
unzip Trimmomatic-0.36.zip

# setup data directory
cd $PROJECT
mkdir data
cd data

# download data from URL
wget https://www.dropbox.com/s/y1ryexxtogx8670/Cace-463_S1_L001_R1_001.fastq.gz
wget https://www.dropbox.com/s/rkd6qr15orcan8i/Cace-463_S1_L001_R2_001.fastq.gz

## check data download: view first few lines of each fastq file
# zcat prints contents of file to screen if file is gzipped
#zcat $R1 | head  
#zcat $R2 | head
