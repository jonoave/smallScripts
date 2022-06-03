#!/usr/bin/env python
# Jun-Hoe, Lee (2021)
# script to generate the md5sum of a file and compare it to the provided md5sum
# designed for sequencing results from DFSZ
## python checkmd5sum.py [list of sequence directory paths]

import os
import sys
import subprocess
import argparse
import configparser

# Import hashlib library (md5 method is part of it)
import hashlib


######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Generates md5sum for R1/R2 files and crosschecks with the md5sum')
parser.add_argument('inputDirList', type=argparse.FileType('r'),
                     help='input file listing all the sample subdirectories containing fastq/')
parser.add_argument('outputCheckmd5sum', type=argparse.FileType('w'),
                     help='output file showing results of cross-checking md5sum')
parser.add_argument('-verbose', action='store_true', default=False,
                    help='turns on verbose mode. Usage: -verbose')  # verbose flag


args = parser.parse_args()

# A. expand list and create symlink to seq files

dirList = args.inputDirlist.readlines()
verboseprint("dirList ", dirList[:5])

filePaths_list = []

for line in dirList:
    R1_path = line + '/fastq/' + line + '_R1.fastqz'
    R2_path = line + '/fastq/' + line + '_R2.fastqz'
    filePaths_list.append(R1_path)
    filePaths_list.append(R2_path)

verboseprint("filePaths_list", filePaths_list[:5])

# B.generate md5suum and crosscheck

for line in filePaths_list:
    7





# File to check
file_name = 'filename.exe'

# Correct original md5 goes here
original_md5 = '5d41402abc4b2a76b9719d911017c592'

# Open,close, read file and calculate MD5 on its contents
with open(file_name) as file_to_check:
    # read contents of the file
    data = file_to_check.read()
    # pipe contents of the file through
    md5_returned = hashlib.md5(data).hexdigest()

# Finally compare original MD5 with freshly calculated
if original_md5 == md5_returned:
    print "MD5 verified."
else:
    print "MD5 verification failed!."#
