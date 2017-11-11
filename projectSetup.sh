#!/bin/bash

## installing software and downloading data

# usage: bash setup.sh URL_R1 URL_R2

# assign variable for location of project directory
PROJECT=`pwd`

# download and install trimmomatic
cd # change to home directory (appropriate for software installation)
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
# if you already have trimmomatic installed, select A to overwrite
unzip Trimmomatic-0.36.zip

# setup data directory
cd $PROJECT
mkdir data
cd data

# download data from URL
wget $1
wget $2

## check data download: view first few lines of each fastq file
# zcat prints contents of file to screen if file is gzipped 
#zcat $R1 | head  
#zcat $R2 | head
# zcat works in GitBash and Linux (e.g., cloud instances), but Macs use gzcat instead
