#!/usr/bin/env python
# -*- coding: latin-1 -*-
# Jun-Hoe, Lee (2021)
# script to rename sample files based on metadata file (tsv) from DataManager (2026)
# input: i.tab-delimited file (txt) converted from an Excel Data Manager downloaed
# ii. location of files to start renaming
# usage:
# python renameFilesDataMana_1.py [list_filenNames] [full_path to directory] -verbose
# WARNING! Make backup before executing script to rename!


import os
import pathlib
import sys
import argparse
import configparser
from tabnanny import verbose
import pandas as pd


######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Renames all files recursively in a directory based on a list of old/new names \
                                         WARNING! Make backup before executing script!')
parser.add_argument('inputFileNames', action="store",
                     help='Measurement file (txt) from Data Manage')
parser.add_argument('dirPath', action='store',
                     help='path to directory of files to rename')
parser.add_argument('-sampleNamesColumn', action='store',
                     help='column name in inputFileNames of the sample names to be replaced',
                     default="Measurement Name")
parser.add_argument('-verbose', action='store_true', default=False,
                    help='turns on verbose mode. Usage: -verbose')  # verbose flag


args = parser.parse_args()

# function to turn on verbose mode
if args.verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print("-v ", arg, ":", args)
else:
    verboseprint = lambda *a: None  # do-nothing function

#########################################################
# A. read in input metadata file

# create dictionary of QBiC_barcodes: given_names
inputFile_df = pd.read_csv(args.inputFileNames, sep='\t', header=0, encoding="utf-16")
verboseprint("inputFile_df ", inputFile_df.head())
             
# Measurement ID* QBiC Sample Id* Sample Name  Sample Pool Group Measurement Name  ...            Flow Cell Sequencing Run Protocol     Index i7     Index i5 Comment
#0  NGSQ2W360001AN-22108713523922     Q2W360001AN         HR2                NaN     26037a001_01  ...  Illumina NovaSeq SP             28+10+10+90  SI-TT-C3-i7  SI-TT-C3-i5     NaN
#1  NGSQ2W360002AW-22109409233281     Q2W360002AW 

## create dictionary of SampleName: QBiC Sample Id*
qbicRename_dict = dict(zip(inputFile_df[args.sampleNamesColumn], inputFile_df['QBiC Sample Id*']))
verboseprint("qbicRename_dict ", list(qbicRename_dict.items())[:5])

# B. get file list 
try:
    # first layer e.g. "AS-630647-LR-56758", "AS-630647-LR-56754"
    file_list = os.listdir(args.dirPath)
    verboseprint("file_list ", file_list)
except:
    print("Incorrect directory path provided. ")
    sys.exit()


for files in file_list:
    for old_name_key in qbicRename_dict.keys():
        if files.find(old_name_key) >= 0: # find substring of old_name_key in file name, if not found it will be -1
            verboseprint("found filename matching sample name", [files, old_name_key])
            # get full path of file, then rename
            oldFilePath = os.path.join(os.path.abspath(args.dirPath), files)
            sample_name = qbicRename_dict.get(old_name_key)
            verboseprint(" sample_name", sample_name)
            # replace only the string of old_name_key in files with sample_name
            newFileName = files.replace(old_name_key, sample_name)
            verboseprint(" newFileName ", newFileName)
            newFilePath = os.path.join(os.path.abspath(args.dirPath), newFileName)
            os.rename(oldFilePath, newFilePath)
        
