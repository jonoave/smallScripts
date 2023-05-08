#!/usr/bin/env python
# Jun-Hoe, Lee (2023)
# script to iterate over each output tsv file of the epitopeprediction/prediction output,
# add in QBIC barcode names and sample condition, filter rows based on certain conditions, 
# input: i. metadata of QBIC barcode, patient group, patient condition
# ii. path to predictions/xxx_prediction_result.tsv
# iii. other conditional arguments e.g. binding, tools etc
# output: tsv file (filtered dataframe with added metadata) 
# usage: 
# python epPred_read_tsv_export.py [input_metadata] [path to predictions folder] [filters] -verbose

import os
import sys
import argparse
import pandas as pd


######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Compile the output from epitope prediction pipeline, \
                            adding sample metadata and other filtering')
parser.add_argument('input_sample_metadata', action="store", type=argparse.FileType('r'),
                     help="tsv: QBiC Code\'\t\'Organism ID\'\t\'Condition: sample_type")               
parser.add_argument('path_abs_to_seq', action='store',
                     help="absolute path to main directory that contains dirs of ep_prediction output files")
parser.add_argument('-suffix', action='store', default = "", 
                     help="Optional uffix to add to each sample with modified tsv, \
                     <QBIC_barcode>_predictions_<suffix>.tsv")
parser.add_argument('-binder', action='store_true', default=False,
                     help="turns on 'binder' column filtering to TRUE only. Usage: -binder")
parser.add_argument('-method', action='store',  nargs='+',
                     help="turns on 'method' (tools) column filtering. Options: netmhcpan-4.1, netmhc-4.0, \
                     and/or syfpeithi-1.0 (separate by space) Default: all 3 methods.")
parser.add_argument('-mergeTables', action='store_true', default=False,
                     help="merge all tsv files into one, using the minimum overlap of columns in all tsv. \
                        Output wil be written as `path_abs_to_seq/merged_predictions.tsv' ")
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
# A. read in input_sample_metadata, store as dictionary

## example of input_sample_metadata,
# QBiC Code	Organism ID	Analyte ID	Condition: sample_type
# QHPRG993AQ	289	FO13192x03_01	blood
# QHPRG994A0	289	FO13192x02_01	primary tumor
# QHPRG995A8	290	FO13191x02_01	primary tumor
# QHPRG996AG	290	FO13191x01_01	blood

# read in readInputSampleMetadata as df
readInputSampleMetadata_df = pd.read_table(args.input_sample_metadata)
# exclude blood samples
readInputSampleMetadata_df = readInputSampleMetadata_df[readInputSampleMetadata_df['Condition: sample_type'] != "blood" ]

# concat columns Organism ID and Condition
readInputSampleMetadata_df['combined'] = readInputSampleMetadata_df['Organism ID'].astype(str) \
+ '_' + readInputSampleMetadata_df['Condition: sample_type']
# create dictionary from both columns
sample_metadata_dict = dict(zip(readInputSampleMetadata_df['QBiC Code'], \
                                readInputSampleMetadata_df['combined']))

verboseprint("sample_metadata_dict   ", list(sample_metadata_dict.items())[:3])


# B. read in predictions_tsv for each sample: 

try:
    dir_list = os.listdir(args.path_abs_to_seq)
    verboseprint("dir_list ", dir_list)
except:
    print("Incorrect directory path provided. ")
    sys.exit()

# walk through each directory and get file names, store in sample_seqFiles_dict
sample_seqFiles_dict = {}
## e.g {"QHPRG00001: ["QHPRG0001_L00x_R1.fastq","QHPRG0001_L00x_R2.fastq"], 
##  "QHPRG00002: ["QHPRG0002_L00x_R1.fastq","QHPRG0002_L00x_R2.fastq"...}

colNames_list = [] # column names placeholder of current merged_pred_df

for dirs in dir_list:
    if not dirs.startswith("VC"):
        continue
    pred_file_name = dirs + "_prediction_result.tsv"
    fullpath_pred_file = os.path.join(args.path_abs_to_seq, dirs, pred_file_name)
    verboseprint("fullpath_pred_file ", pred_file_name)

    # read in as tsv with pandas 
    ep_pred_df = pd.read_table(fullpath_pred_file, dtype = {'chr': str})
    verboseprint("dimension of just read ep_pred_df ", ep_pred_df.shape)
    verboseprint(" ep_pred_df head ",  ep_pred_df.head())

    # append metadata of QBIC barcode, patient_group, and sample_status
    ## get file name:
    qbic_barcode = dirs.split("_")[0][4:] # starts with VC01
    verboseprint("qbic_barcode from directory: ", qbic_barcode)
    patient_group, sample_status = sample_metadata_dict[qbic_barcode].split("_")

    ep_pred_df['QBIC_barcode'] = qbic_barcode
    ep_pred_df['patient'] = patient_group
    ep_pred_df['condition'] = sample_status

    # move columns to the start of the df
    cols_to_move = ['QBIC_barcode', 'patient', 'condition']
    ep_pred_df = ep_pred_df[cols_to_move + [ col for col in ep_pred_df.columns if col not in cols_to_move ]]
    verboseprint("added qbic barcode to ep_pred_df ", ep_pred_df.head())
    
    # apply filterings
    if args.binder:
        ep_pred_df = ep_pred_df[ep_pred_df["binder"] == True]

    if args.method:
        try:
            for tool in args.method:
                ep_pred_df= ep_pred_df[ep_pred_df["method"] == tool]
        except:
            print("-method argument(s) provided do not exist in the predicted tsv output")
            sys.exit()
    
    # write modifed df as <qbic_barcode>_modified_ts
    ## replace NAN values with empty strings
    ep_pred_df = ep_pred_df.fillna("")
    # sort by columns "QBIC_barcode", then "chr", "length"
    ep_pred_df  = ep_pred_df.sort_values(["QBIC_barcode", "chr", "length"])
    verboseprint("dimension of ep_pred_df to be written ", ep_pred_df.shape)
    outfile_ep_pred_df = qbic_barcode +  "_prediction_result_" + args.suffix + ".tsv"
    fullpath_outfile_ep_pred_df  = os.path.join(args.path_abs_to_seq, dirs, outfile_ep_pred_df)
    ep_pred_df.to_csv(fullpath_outfile_ep_pred_df, sep="\t", index=False) 

    # create a common merged tables if args.mergeTables
    if args.mergeTables:
        if not colNames_list: # for first sample
            colNames_list = ep_pred_df.columns.tolist()
            verboseprint("colnames for first merged_ep_pred_df ", colNames_list)
            merged_ep_pred_df = ep_pred_df
            # keepc current ep_pred_df as is
        else: # check colnames against previous samples (colNames_list)
            # current_colNames_list = ep_pred_df.columns
            cols_ep_pred_df_set = set(ep_pred_df.columns.tolist())
            verboseprint("Dimensions of ep_pred_df ", ep_pred_df.shape)
            cols_merged_ep_pred_df_set = set(merged_ep_pred_df.columns.tolist())
            verboseprint("Dimensions of merged_ep_pred_df, ", merged_ep_pred_df.shape)
            overlap_colNames_set = cols_ep_pred_df_set.intersection(cols_merged_ep_pred_df_set)
            # make sure both current merged_ep_pred_df and the new ep_pred_df have the same columns as the overlap
            # remove extra columns that are not found in overlap
            if len(cols_ep_pred_df_set) > len(overlap_colNames_set): # more columns in ep_pred_df
                complement_ep_pred_df_set = cols_ep_pred_df_set.difference(overlap_colNames_set)
                verboseprint(" colNames ep_pred_df before dropping columns ", ep_pred_df.columns.tolist()[:7])
                verboseprint("complement_ep_pred_df_set ", complement_ep_pred_df_set )
                mod_ep_pred_df = ep_pred_df[ep_pred_df.columns.difference(list(complement_ep_pred_df_set), sort = False)]
                verboseprint(" colNames mod_ep_pred_df after dropping columns ", mod_ep_pred_df.columns.tolist()[:7])
            if len(cols_merged_ep_pred_df_set) > len(overlap_colNames_set): # more columns in merged_ep_pred_df
                complement_merged_ep_pred_df_set = cols_merged_ep_pred_df_set.difference(overlap_colNames_set)
                mod_merged_ep_pred_df = merged_ep_pred_df[merged_ep_pred_df.columns.difference(list(complement_merged_ep_pred_df_set), sort = False)]
            else:
                # the columns in overlap_colNames_set are the same with merged_ep_pred_df and ep_pred_df
                mod_merged_ep_pred_df = merged_ep_pred_df.copy()
                mod_ep_pred_df = ep_pred_df
            
            # append mod_ep_pred_df to mod_merged_ep_pred_df
            # just to be safe, make sure the order of columns in ep_pred_df = merged_ep_df           
            mod_merged_ep_pred_df = mod_merged_ep_pred_df[mod_ep_pred_df.columns.tolist()]
            # mod_ep_pred_df = mod_ep_pred_df[mod_merged_ep_pred_df.columns.tolist()]
            merged_ep_pred_df = pd.concat([mod_merged_ep_pred_df, mod_ep_pred_df])
            verboseprint('Dimenions of concatenated merged_ep_pred_df ', merged_ep_pred_df.shape)
            # reset colNames_list 
            colNames_list = merged_ep_pred_df.columns.tolist()

# write out mergeTables 
if args.mergeTables:
    merged_ep_pred_df = merged_ep_pred_df.fillna("")
    # sort by column "QBIC_barcode", then "chr", "length"
    merged_ep_pred_df  = merged_ep_pred_df.sort_values(["QBIC_barcode", "chr", "length"])
    verboseprint("dimension of merged_ep_pred_df (concat) to be written \n ", merged_ep_pred_df.shape)
    fullpath_outfile_merged_ep_pred_df = os.path.join(args.path_abs_to_seq, "merged_ep_pred_df.tsv")

    merged_ep_pred_df.to_csv(fullpath_outfile_merged_ep_pred_df, sep="\t", index=False)

    



                                                                       
                                                                       

    






