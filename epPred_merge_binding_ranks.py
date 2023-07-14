#!/usr/bin/env python
# Jun-Hoe, Lee (2023)
# child script to epPred_read_tsv_export.py
# goes through the <..>_binding_rank.tsv file of each sample, and concatenate them into one table. 
# the first five columns are common ["QBIC_barcode","patient","condition","sequence", "gene", ],
# but the subsequent 6 columns of HLAs vary by sample
# input: main directory, where the subdirectories contain a <..>binding_rank.tsv file
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

parser = errorDisplayParser(description='After running epPred_read_tsv_export.py to generate <..>_binding_rank.tsv \
                            for each sample, now combine them into one merged table')
parser.add_argument('path_abs_to_seq', action='store',
                     help="absolute path to main directory that contains dirs of bindRankiction output files. \
                        Subdirectories also contains the <..>_binding_rank.tsv file")
parser.add_argument('dirs_match', action='store',
                     help="QBIC barcode prefix (first 5 characters), e.g. QHPRG.")
# parser.add_argument('output_file', action="store", type=argparse.FileType('w'),
#                      help='output file name of the merged <..>_binding_rank.tsv files ')
parser.add_argument('-verbose', action='store_true', default=False,
                    help="turns on verbose mode. Usage: -verbose")  # verbose flag


args = parser.parse_args()

# function to turn on verbose mode
if args.verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print("-v ", arg, ":", args)
else:
    verboseprint = lambda *a: None  # do-nothing function

#########################################################

# A. read in <..>_binding_rank.tsv for each sample: 

try:
    dir_list = os.listdir(args.path_abs_to_seq)
    verboseprint("dir_list ", dir_list)
except:
    print("Incorrect directory path provided. ")
    sys.exit()

# filter only for directories containing QBIC_barcode
sample_dir_list = [i for i in dir_list if args.dirs_match in i]


commonCols_list = [] # column names placeholder of current merged_pred_df

for dirs in sample_dir_list:
    qbic_match = args.dirs_match + "\d{3}\w{2}"
    dirs_search = re.search(qbic_match, dirs)
    if not dirs_search:
        print("Cannot find QBIC sample matching project code: ", dirs)
        continue
    # get QBIC barcode of sample, sometimes it starts with VCxxQHPRGxxxxx
    qbic_sample_portion  = dirs.split("_")[0]
    qbic_find_sample = re.search(qbic_match, qbic_sample_portion)
    qbic_sample_ID = qbic_find_sample.group()
    verboseprint("qbic_sample_ID  ", qbic_sample_ID)
    bindingRank_file_name = qbic_sample_ID + "_binding_rank_.tsv"
    fullpath_bindingRank_file = os.path.join(args.path_abs_to_seq, dirs, bindingRank_file_name)
    verboseprint("fullpath_bindingRank_file :  ", fullpath_bindingRank_file)

    # read in as tsv with pandas
    try:
        bindRank_df = pd.read_table(fullpath_bindingRank_file, dtype = {'chr': str})
        verboseprint("dimension of just read bindRank_df ", bindRank_df.shape)
        verboseprint(" bindRank_df head ",  bindRank_df.head())
    except:
        print("No binding_rank.tsv file found for:  ", bindingRank_file_name)
        continue
    
    # create dynamic updating tables (for each sample) for merge_tables and HLA_distribution
    # A. Compiling information for common merged tables if args.mergeTables
    if not commonCols_list: # first table, where it is empty
        commonCols_list = bindRank_df.columns.tolist()[:7] # only the first 7 columns are common
        verboseprint("common cols for merged table ", commonCols_list)
        merged_bindRank_df = bindRank_df
        # keepc current bindRank_df as is
    else: # merge columns based on commonCols_list
        verboseprint('Dimenions of bindRank_df before joining ', bindRank_df.shape)
        verboseprint('Dimenions of merged_bindRank_df before joining ', merged_bindRank_df.shape)
        mod_merged_bindRank_df = merged_bindRank_df.copy()
        merged_bindRank_df = pd.concat([mod_merged_bindRank_df, bindRank_df])
        verboseprint('Dimenions of merged_bindRank_df after joining ', merged_bindRank_df.shape)
        verboseprint("merged_bindRank_df head ", merged_bindRank_df.head())

# B. sort columns before writing them out.
colsNames_list = merged_bindRank_df.columns.tolist()
verboseprint("colsNames_list ", colsNames_list)
# check that the commonCols_list is still found within colsNames_list
colsNames_set = set(colsNames_list)
verboseprint("colsNames_set ", colsNames_set )
commonCols_set = set(commonCols_list)
verboseprint("commonCols_set ", commonCols_set )

if (commonCols_set.issubset(colsNames_set)):
    print("yes, it is a set")
    hla_set = commonCols_set.difference(colsNames_set)
    hla_set = colsNames_set.difference(commonCols_set)
    verboseprint("hla_set" , hla_set)
    hla_set_2 = colsNames_set - commonCols_set
    verboseprint(hla_set_2, hla_set_2)
    # sort hla_set alphabetically
    hla_list = sorted(list(hla_set))
    verboseprint("hla_list, sorted ", hla_list)
else:
    print("Common first five columns ['QBIC_barcode','patient','condition','sequence', 'gene'] \
          not found in final merged_binding_rank df ")
    sys.exit()

# reorder columns in merged_bindRank_df 
colsNames_sorted_list = commonCols_list + hla_list
merged_bindRank_sorted_df = merged_bindRank_df[colsNames_sorted_list]
verboseprint('Dimensions of merged_bindRank_sorted_df after reordering ', merged_bindRank_sorted_df.shape)

# write out final table
fullpath_outfile_merged_bindRank_df = os.path.join(args.path_abs_to_seq, "merged_binding_rank_HLA.tsv")
merged_bindRank_sorted_df.to_csv(fullpath_outfile_merged_bindRank_df, sep="\t", na_rep = "NaN", index=False)        


