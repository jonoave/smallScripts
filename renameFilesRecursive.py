#!/usr/bin/env python
# Jun-Hoe, Lee (2021)
# script to rename all files recursively in a folder based on input list
# input: i.tsv of old names and new names
# ii. location to start renaming
# usage:
# python renameFilesRecursive.py [list_filenNames] [working directory] -verbose
# WARNING! Make backup before executing script to rename!


import os
import pathlib
import sys
import subprocess
import argparse
import configparser


######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Renames all files recursively in a directory based on a list of old/new names \
                                         WARNING! Make backup before executing script!')
parser.add_argument('inputFileNames', action="store", type=argparse.FileType('r'),
                     help='tsv of old and new file names')
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

# create dictionary of oldnames: new names:
readInputFileNames = args.inputFileNames.readlines()
qbicRename_dict = {}
for line in readInputFileNames:
    line = line.strip() # remove \n
    splitLine = line.split("\t")
    qbicRename_dict[splitLine[1]] = splitLine[0]

verboseprint("qbicRename_dict ", list(qbicRename_dict.items())[:5])

# B. rename files and folder recursively in directory
## do this in 2 steps to avoid breaking
## first rename all files
try:
    # first layer e.g. "AS-630647-LR-56758", "AS-630647-LR-56754"
    dir_list = os.listdir(args.dirPath)
    verboseprint("dir_list ", dir_list)
except:
    print("Incorrect directory path provided. ")
    sys.exit()


for dirs in dir_list:
    if dirs in qbicRename_dict:
        verboseprint(" dirs in qbicRename_dict ", dirs)
        currentPath = pathlib.Path().resolve()
        fullpath = os.path.join(currentPath , dirs)
        verboseprint("fullpath ", fullpath)
        qbicBarcode = qbicRename_dict.get(dirs)
        for subdirs, subSubdirs, files in os.walk(fullpath):
            verboseprint("subdirs in os.walk(dirs)" , subdirs)
            # second layer, e.g. fastq/xxx_R1.fastqz
            #qbicBarcode = qbicRename_dict.get(dirs)
            for fileName in files:
                verboseprint('fileName', fileName)
                if fileName.find(dirs) >= 0: # find substring of dirs in filename, if not found it will be -1
                    print('found dirs in filename')
                    subdirectoryPath = os.path.relpath(subdirs, dirs) #get the path to your subdirectory
                    verboseprint("subdirectoryPath ", subdirectoryPath)
#                    filePath = os.path.join(subdirectoryPath, fileName) #get the path to your file
                    oldFilePath = os.path.join(fullpath, subdirectoryPath, fileName)
                    verboseprint("oldFilePath",oldFilePath)
                    newFileName = fileName.replace(dirs, qbicBarcode) #create the new name
                    newFilePath = os.path.join(fullpath, subdirectoryPath, newFileName)
                    # get path of old and new file names
                    os.rename(oldFilePath, newFilePath) #rename the file

        # now rename the subdirs
        if os.path.basename(fullpath) == dirs:
            newFullPath = fullpath.replace(dirs, qbicBarcode)
            os.rename(fullpath, newFullPath)
