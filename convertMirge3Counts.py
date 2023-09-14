#!/usr/bin/env python
# Jun-Hoe, Lee (2023)
# script to read in miRNA counts mir.Counts.csv output file from mirge3 workflow , 
# It's a simple csv file that needs to be converted to tsv, the gene (mirna) columns duplicated,
# and the headers renamed
# input: <miR.Counts.csv) 
# output: <..>_merged_binding_rank. tsv file (filtered dataframe with added metadata) 
# usage: 
# python epPred_read_tsv_export.py [input_metadata] [path to predictions folder] [filters] -verbose

import os
import sys
import argparse
import pandas as pd
import re



######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()   
        sys.exit(2)

parser = errorDisplayParser(description='Convert the output of nf-core/smrnaseq pipeline in edgeR, ie. the <hairpin,mature> \
                            count.tsv files into a merged_gene_count.tsv file as input for rnadeseq')
parser.add_argument('input_counts_csv', action='store',
                    help="miR.Counts.csv file from mirge3 output")
parser.add_argument('-output_file', action="store", type=argparse.FileType('w'),
                    default='merged_gene_counts.tsv', help="output file name of processed merged mirna counts. \
                    Default. 'merged_gene_counts.tsv'.")
parser.add_argument('-verbose', action='store_true', default=False,
                    help="turns on verbose mode. Usage: -verbose")  # verbose flag

args = parser.parse_args()

# function to turn on verbose mode
if args.verbose:
    def verboseprint(*args, **kwargs):
        print("-v ", *args, **kwargs)
else:
    verboseprint = lambda *a, **k: None # do-nothing function
######################################################
# A.read in input files:

# read in, ignoring rownames
mirnaCount_df = pd.read_csv(args.input_counts_csv)
verboseprint(mirnaCount_df.head())
column_names = list(mirnaCount_df.columns)
verboseprint("column names \n", column_names)

## convert decimals to integer
keep_columns = mirnaCount_df.columns.drop('miRNA') # exclude first column
# verboseprint("keep_columns ", keep_columns)
mirnaCount_df[keep_columns] = mirnaCount_df[keep_columns].round(0).astype(int)
verboseprint("convert decimals to integer \n",mirnaCount_df.head())


# duplicate first column "miRNA" as "gene_name" 
mirnaCount_df['gene_name'] = mirnaCount_df.loc[:, 'miRNA']
verboseprint("duplicate first column \n", mirnaCount_df.head())

# rename column 1 from "Geneid"
mirnaCount_df = mirnaCount_df.rename(columns={"miRNA": "Geneid"})

# move column "gene_name" from last to second position
genename_column = mirnaCount_df.pop('gene_name')
mirnaCount_df.insert(1, 'gene_name', genename_column)
verboseprint("mirnaCount_df, moved column genename  \n", mirnaCount_df.head())

# write the output of the df 
mirnaCount_df.to_csv(args.output_file, sep="\t", index=False )


