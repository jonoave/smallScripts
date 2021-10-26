#!/usr/bin/env python
# Jun-Hoe, Lee (2020)
# script to parse the output of clean/filtered output of changeO /input of igphyml
# to generate a fasta alignment file 
# usage:
## python getIgphymlMSAclone.py [clone_id] [xxx_db-pass_productive-T_clone-pass_germ-pass.tsv] [output MFA file] -verbose


import os
import sys
import subprocess
import argparse
import configparser

# from Bio import SeqIO
import pandas as pd

######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Generates a multiple sequence aligment including germline sequence from a clone cluster/group')
parser.add_argument('cloneID', action="store", type=int,
                     help='cloneID number')
parser.add_argument('inputDB_pass', action='store',
                     help='input file from changeO pipeline, \
                     \ e.g. xxx_db-pass_productive-T_clone-pass_germ-pass.tsv')
parser.add_argument('outputMSA', type=argparse.FileType('w'),
                     help='output file to write the generated multiple sequnce alignment')
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

# A. read in input files
# read in as panda
try:
    germSeq_ID = str(args.cloneID) + "_GERM"
    df = pd.read_csv(args.inputDB_pass, sep='\t')
    verboseprint("df head", df.head())
except:
    print("Cannot read input file from changeO pipeline")

## subset df based on clone_id
df_clone = df.loc[df['clone_id'] == args.cloneID]
verboseprint("df_clone head ", df_clone.head())

# get germline sequence
## grab the germline_alignment_d_mask from the first clone_id

germSeq = df_clone["germline_alignment_d_mask"].iloc[0]
germSeq_write = ">" + germSeq_ID + '\n' + germSeq
verboseprint("germSeq_write ", germSeq_write)
args.outputMSA.write(germSeq_write)


# create lists from sequence name  ('sequence_id') and sequence ('germline_alignment_d_mask')
seqName_list = df_clone["sequence_id"].tolist()
df_clone = df_clone.copy(deep=False)  # Ensuring a copy is made
df_clone['sequence_alignment'] = df_clone['sequence_alignment'].str.replace('.','-')
sequence_list = df_clone["sequence_alignment"].tolist()

for seqName, seq in zip(seqName_list, sequence_list):
    writeLine = '\n' + ">" + seqName + '\n' + seq
    args.outputMSA.write(writeLine)










# combine seqID and sequence
