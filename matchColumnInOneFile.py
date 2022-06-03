#!/usr/bin/env python
# Jun-Hoe, Lee (2022)
# script to read in file1, file2. And sort rows in file2 based on file1`` 
# key to sort file1, file 2 should be column 1
# input: i. file1 (tsv) only 1 column (desired sort)
# ii. file2 (tsv) multiple columns allowed
# usage:
# python matchColumnInOneFile.py [file1] [file2] [outputfile] -verbose

import os
import pathlib
import sys
import subprocess
import argparse
import configparser

import pandas as pd


######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Match row order in file2 to those of file1, based on col1 in both files')
parser.add_argument("input_file1", action="store", type=argparse.FileType('r'),
                     help='file1 containing 1 column (desired sort')
parser.add_argument("input_file2", action="store", type=argparse.FileType('r'),
                     help='file2 containing multiple columns, column1 contents (unsorted) match file1')
parser.add_argument("output_file", action="store", type=argparse.FileType('w'),
                     help='output file ')
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
# A. read in input_file1

df_file1 = pd.read_csv(args.input_file1, sep='\t', header=None)
verboseprint("head df_file1 ", df_file1.head())
# convert column to list
file1_list = df_file1.iloc[:, 0].tolist()
verboseprint("file1_list ", file1_list[:5])

# B. read in input_file2
df_file2 = pd.read_csv(args.input_file2, sep='\t', header=None)
df_file2.columns = ["sorter_key", "QBIC barcode", "Secondary name"]
verboseprint("head df_file2 ", df_file2.head())

# C. sort rows in df_file2 based on file1_list
 ## Create a dummy df with the required list and the col name to sort on
dummy = pd.Series(file1_list, name = "sorter_key").to_frame()
 ## Use left merge on the dummy to return a sorted df
df_sorted = pd.merge(dummy, df_file2, on = 'sorter_key', how = 'left')
verboseprint("head df_sorted ", df_sorted.head())

## sort based on column in file_1

# D. write output
df_sorted.to_csv(args.output_file, sep="\t")




