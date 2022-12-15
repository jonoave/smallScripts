#!/usr/bin/env python
# Jun-Hoe, Lee (2020)
# script to parse the output of clean/filtered output of changeO /input of igphyml
# to generate a fasta alignment file 
# usage:
## python getIgphymlMSAclone.py [xxx_db-pass_productive-T_clone-pass_germ-pass.tsv] [xxx_db-pass_productive-T_clone-pass_germ-pass.tsv] [output MFA file] -verbose


import os
import sys
import subprocess
import argparse
import configparser

# from Bio import SeqIO
from ete3 import Tree
import pandas as pd

######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Generates a multiple sequence aligment including germline sequence from a clone cluster/group')
parser.add_argument('inputDB_pass_tsv', action='store',
                     help='input file from changeO pipeline, \
                     e.g. xxx_db-pass_productive-T_clone-pass_germ-pass.tsv')
parser.add_argument('inputDB_phylo_tree', action='store',
                     help='input of newick tree (transformed output from BuildTree), \
                     e.g. xxx_phyloTree_newick.txt')
parser.add_argument('outputMSA', type=argparse.FileType('w'),
                     help='output fasta file to write the generated multiple sequnce alignment')
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

# A. read in newick tree, get clone sequence names
verboseprint("args.inputDB_phylo_tree", args.inputDB_phylo_tree)
try:
    clone_tree = Tree(args.inputDB_phylo_tree)
except:
    print("Cannot read input newick tree")
    sys.exit()

# get list of leaves (clone sequence names)
clone_seqs_list = clone_tree.get_leaf_names()
verboseprint("clone_seqs_list", clone_seqs_list[:5])
verboseprint("clone_seqs_list", clone_seqs_list[-5])
verboseprint("length clone_seqs_list ", len(clone_seqs_list))
# in the inputDB_pass_tsv file, seq names are originally with ":"
clone_seq_repl_list = [seqs.replace("-",":") for seqs in clone_seqs_list]
verboseprint("clone_seq_repl_list", clone_seq_repl_list[:5])


# B. read in input files
# read in xx_germ-pass.tsv into pandas
try:
    df = pd.read_csv(args.inputDB_pass_tsv, sep='\t')
    verboseprint("df head", df.head())
except:
    print("Cannot read input file from changeO pipeline")
    sys.exit()
    
# subset df based on clone_seq_repl_list, i.e. only those used in tree reconstruction
clone_seqs_df = df[df['sequence_id'].isin(clone_seq_repl_list)]
## check number of rows
verboseprint("number of rows in clone_seqs_df ", clone_seqs_df.shape[0]) 

# get germline sequence 
## grab the germline_alignment_d_mask from the first clone sequence
germ_seq = clone_seqs_df.iloc[1]['germline_alignment_d_mask']
verboseprint("germ_seq  ", germ_seq)
## grab the clone_ID of the germline being used
germ_seq_clone_ID = clone_seqs_df.iloc[1]['clone_id']
verboseprint("germ_seq_clone_ID ", germ_seq_clone_ID)

# create lists from sequence name  ('sequence_id') and sequence ('sequence_alignment'),
# because we want to maintain the order of `sequence_id'. paired with their sequences
clone_seq_name_list = clone_seqs_df["sequence_id"].tolist()
clone_seq_list = clone_seqs_df['sequence_alignment'].tolist()
## replace dots in sequences with "-" as gaps
clone_seq_repl_list = [seqs.replace(".","-") for seqs in clone_seq_list]


# C. write output .fas file
# write germline sequence first
germ_seq_name = str(germ_seq_clone_ID) + "_GERM"
verboseprint("germ_seq_name ", germ_seq_name)
writeLine_germseq = ">" + germ_seq_name + '\n' + germ_seq 
args.outputMSA.write(writeLine_germseq)

# write the remainin sequences.
for seqName, seq in zip(clone_seq_name_list, clone_seq_repl_list):
    writeLine = '\n' + ">" + seqName + '\n' + seq
    args.outputMSA.write(writeLine)
