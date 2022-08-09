#!/usr/bin/env python
# Jun-Hoe, Lee (2021)
# script to rename all files recursively in a folder based on input list
# folders are already QBiC barcodes-named, but the files within aren't
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
    qbicRename_dict[splitLine[0]] = splitLine[1]

verboseprint("qbicRename_dict ", list(qbicRename_dict.items())[:5])

# B. rename files within each folder
try:
    # first layer e.g. "AS-630647-LR-56758", "AS-630647-LR-56754"
    dir_list = os.listdir(args.dirPath)
    verboseprint("dir_list ", dir_list)
except:
    print("Incorrect directory path provided. ")
    sys.exit()


for dirs in dir_list:
    if dirs in qbicRename_dict.keys(): # directory matching QBiC barcode
        verboseprint(" dir with QbiC barcode in qbicRename_dict ", dirs)
        # get full path of directory
        dir_fullpath = os.path.join(os.path.abspath(args.dirPath), dirs)
        sample_name = qbicRename_dict.get(dirs)
        verboseprint(" sample name: ", sample_name) 
        # walk directory and get the names of files inside
        for files in os.listdir(dir_fullpath):
            if files.find(sample_name) >= 0: # find substring of sample_name in file name, if not found it will be -1
                verboseprint("found filename matching sample name", files)
                # get full path of file, then rename
                oldFilePath = os.path.join(dir_fullpath, files)
                newFileName = files.replace(sample_name, dirs)
                verboseprint( " newFileName ", newFileName)
                newFilePath = os.path.join(dir_fullpath, newFileName)
                os.rename(oldFilePath, newFilePath)