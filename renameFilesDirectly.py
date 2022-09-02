#!/usr/bin/env python
# Jun-Hoe, Lee (2021)
# script to rename all files directly, in a folder based on input list
# input: i.tsv of new (QBIC) and old names (e.g. DKFZ names)
# ii. location to start renaming
# usage:
# python renameFilesUnderQBiC.py [list_filenNames] [full_path to directory] -verbose
# WARNING! Make backup before executing script to rename!


import os
import pathlib
import sys
import argparse
import configparser
from tabnanny import verbose


######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Renames all files recursively in a directory based on a list of old/new names \
                                         WARNING! Make backup before executing script!')
parser.add_argument('inputFileNames', action="store", type=argparse.FileType('r'),
                     help='tsv file: QBIC_barcodes\'\t\'Old_names')
parser.add_argument('dirPath', action='store',
                     help='path to directory')
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
# A. read in input file Renames

# create dictionary of QBiC_barcodes: given_names
## the QBiC_barcodes are also the folder names 
readInputFileNames = args.inputFileNames.readlines()
qbicRename_dict = {}
for line in readInputFileNames:
    line = line.strip() # remove \n
    splitLine = line.split("\t")
    qbicRename_dict[splitLine[1]] = splitLine[0]

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
        
