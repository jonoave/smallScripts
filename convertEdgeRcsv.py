#!/usr/bin/env python
# Jun-Hoe, Lee (2023)
# script to read in miRNA counts <hairpin/mature>_counts.csv file from nf-core/smrnaseq , 
# where the rownames are the sample names. Transpose it and convert it to input for rnadeseq
# input:  
# output: merged_gene_counts.tsv  
# usage: 
# python convertEdgeRcsv.py [<hairpin/mature>_counts.csv] [-options]

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
parser.add_argument('input_counts_csv', nargs='+', action='store',
                    help="<hairpin/mature>_counts.csv files, separated by space")
parser.add_argument('-output_file', action="store", type=argparse.FileType('w'),
                    default='merged_gene_counts.tsv', help="output file name of tranposed/processed/merged mirna counts. \
                    Default. 'merged_gene_counts.tsv'.")
parser.add_argument('-old_v1_pipeline', action='store_true', default=False,
                    help="old pipeline version, which has mirna tags are now used in place of gene_names")  
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
# function to transpose dataframe:
def transpose_df(input_df):
    transposed_df = input_df.transpose() 
    # the transposed df has a new header index; remove this.
    transposed_df = transposed_df.rename(columns=transposed_df.iloc[0]).drop(transposed_df.index[0])
    # re-order columns alphabetically 
    transposed_df = transposed_df.reindex(sorted(transposed_df.columns), axis=1)
    verboseprint(transposed_df.head())
    # show number of rows (mrnas)
    verboseprint("number of row in transposed df: ", transposed_df.shape[0])


    return transposed_df

######################################################

# A.read in input files:

mirnaType_list = [] # list to store mirnatypes (hairpin/mature)
df_list = [] # list to store the read in dfs
df_samples_list = [] # list to store values in column 1 (sample_names)


for file_name in args.input_counts_csv:
    verboseprint(file_name)
    # get first name (hairpin_or mature)
    mirnaType_list.append(file_name.split("_")[0])
    verboseprint(mirnaType_list)

    # read in as csv with pandas
    mirnaCount_df = pd.read_csv(file_name)
    verboseprint(mirnaCount_df.head())
  
    # remove mirna suffix from sample nmames in first column, e.g QSCNN002A5.mature , QSCNN002A5.hairpin
    mirnaCount_df.iloc[:, 0] =  mirnaCount_df.iloc[:, 0].apply(lambda x: x.split(".")[0])
    # now get sample names, for subsequent check that the same samples are present in all dfs
    samples_mirna_list = mirnaCount_df.iloc[:, 0].tolist()
    verboseprint(samples_mirna_list)
    df_samples_list.append(samples_mirna_list)
    verboseprint(df_samples_list)

    # store df in list
    df_list.append(mirnaCount_df)

tranposed_df_list = [] # list to store tranposed df
if len(args.input_counts_csv) > 1: # more than one input files provided    
    # check whether the row names (sample names) are the same in both/all files
    for samples in df_samples_list[1:]:
        if set(samples) != set(df_samples_list[0]):
            print("mirna count files have different sample names - check rows!")
            sys.exit()

    # transpose both dfs 
    for dfs in df_list:
        transposed_df = transpose_df(dfs)
        tranposed_df_list.append(transposed_df)
    
    #  merge df
    merged_df = pd.concat(tranposed_df_list)

else: # just 1 sample, so just transpose df
    merged_df = transpose_df(df_list[0])

# modify the merged_df
## i. split first column (actually index) values into "Geneid" and "gene_name"
## e.g "hsa-let-7a-1_MI0000060_Homo_sapiens_let-7a-1_st.." take the first 2 values split by "_"

# reset index, change into a column
merged_df.index.name = 'Geneid'
merged_df = merged_df.reset_index()
verboseprint("merged_df renamed index to Geneid \n", merged_df.head())
# duplicate column "Geneid"
merged_df['gene_name'] = merged_df.loc[:, 'Geneid']
verboseprint("merged_df, create column genename  \n", merged_df.head())
# move column "gene_name" from last to second position
genename_column = merged_df.pop('gene_name')
merged_df.insert(1, 'gene_name', genename_column)
verboseprint("merged_df, moved column genename  \n", merged_df.head())

## the new pipeline v > 2.0 only shows the mirID, no mirTag. 
## so we skip the following step. 
##  now edit the values in columns 'Geneid' and 'gene_name'
if args.old_v1_pipeline:
    merged_df['Geneid'] =  merged_df['Geneid'].apply(lambda x: x.split("_")[0])
    merged_df['gene_name'] =  merged_df['gene_name'].apply(lambda x: x.split("_")[1])
    verboseprint("edited values in Geneid and genename  \n", merged_df.head())

# write the output of the df 
merged_df.to_csv(args.output_file, sep="\t", index=False )
