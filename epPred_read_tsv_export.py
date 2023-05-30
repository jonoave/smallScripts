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
import re

from collections import OrderedDict, defaultdict


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
parser.add_argument('dirs_match', action='store',
                     help="QBIC barcode prefix (first 5 characters), e.g. QHPRG.")
parser.add_argument('-suffix', action='store', default = "", 
                     help="Optional suffix to add to each sample with modified tsv, \
                     <QBIC_barcode>_predictions_<suffix>.tsv")
parser.add_argument('-binder', action='store_true', default=False,
                     help="turns on 'binder' column filtering to TRUE only. Usage: -binder")
parser.add_argument('-method', action='store',  nargs='+',
                     help="turns on 'method' (tools) column filtering. Options: netmhcpan-4.1, netmhc-4.0, \
                     and/or syfpeithi-1.0 (separate by space) Default: all 3 methods.")
parser.add_argument('-mergeTables', action='store_true', default=False,
                     help="merge all tsv files into one, using the minimum overlap of columns in all tsv. \
                        Output wil be written as `path_abs_to_seq/merged_predictions.tsv'. This also creates \
                        2 tsv files that compiles number of HLAs by SB and WB, respectively by 'primary tumour' \
                        and 'metastasis'. ")
parser.add_argument('-longest', action='store_true', default=False,
                     help="only keep the longest predicted peptide sequence for each chr/position.")
parser.add_argument('-bindStrength', action='store_true', default=False,
                     help="compiles epitopes from all samples that are ranked as SB or WB, by looking at the \
                          'HLAxx:xx rank' column, output as 'epitope_bindStrength.tsv`. Usage: -bindStrength")
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
# function to update counter (value) of a key in dict
def update_dict_counter(counter_dict, counter_dict_key):
    if counter_dict_key in counter_dict:
        counter_dict[counter_dict_key] += 1
    else:
        counter_dict[counter_dict_key] = 1

def update_global_hla_counter(counter_dict, counter_dict_key, binding_str, binding_num):
    verboseprint(" counter_dict_key, binding ",  [counter_dict_key, binding_str])
    verboseprint(" counter_dict before updating",  list(counter_dict.items()))
    if counter_dict_key not in counter_dict:
        if binding_str== "SB":
            counter_dict[counter_dict_key] = [binding_num, 0]
        else:
            counter_dict[counter_dict_key] = [0, binding_num]
    else:
        binding_counter_list = counter_dict.get(counter_dict_key)
        if binding_str == "SB":
            binding_counter_list[0] += binding_num
        elif binding_str== "WB":
            binding_counter_list[1] += binding_num
        verboseprint("updated binding_counter_list ", binding_counter_list)
    # update dictionary with binding_counter value
    # updated_dict = {counter_dict_key: binding_counter_list}
    # counter_dict.update(updated_dict)

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

# dictionaries to store binding strength of epitopes by HLA across samples
ep_binding_hlas_dict = {}
## store as "condition"_"WB/SB"_"HLAxxx"": count

for dirs in dir_list:
    # qbic_match = re.escape(args.dirs_match + "\d{3}\w\d")
    qbic_match = args.dirs_match + "\d{3}\w{2}"

    verboseprint("qbic_match ", qbic_match)
    dirs_search = re.search(qbic_match, dirs)
    if not dirs_search:
        print("Cannot find predictions file matching: ", dirs)
        continue
    verboseprint("dirs_search ", dirs_search.group())
    pred_file_name = dirs + "_prediction_result.tsv"
    fullpath_pred_file = os.path.join(args.path_abs_to_seq, dirs, pred_file_name)
    verboseprint("fullpath_pred_file ", pred_file_name)

    # read in as tsv with pandas
    try:
        ep_pred_df = pd.read_table(fullpath_pred_file, dtype = {'chr': str})
        verboseprint("dimension of just read ep_pred_df ", ep_pred_df.shape)
        verboseprint(" ep_pred_df head ",  ep_pred_df.head())
    except:
        print("No prediction_result.tsv file found for:  ", dirs)
        continue

    # append metadata of QBIC barcode, patient_group, and sample_status
    ## get file name:
    # qbic_barcode = dirs.split("_")[0][4:] # starts with VC01
    qbic_barcode = dirs_search.group() # starts with VC01
    verboseprint("qbic_barcode from directory: ", qbic_barcode)
    patient_group, sample_status = sample_metadata_dict[qbic_barcode].split("_")

    ## create following columns that correspondes to length of ep_pred_df
    ep_pred_df['QBIC_barcode'] = qbic_barcode
    ep_pred_df['patient'] = patient_group
    ep_pred_df['condition'] = sample_status

    # move columns to the start of the df
    cols_to_move = ['QBIC_barcode', 'patient', 'condition']
    ep_pred_df = ep_pred_df[cols_to_move + [ col for col in ep_pred_df.columns if col not in cols_to_move ]]
    verboseprint("added qbic barcode to ep_pred_df ", ep_pred_df.head())

    # drop columns CSQ_IMGAG to save space
    ep_pred_df = ep_pred_df.drop(columns= ["CSQ_IMGAG"]) 
    
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

    # parameter to keep only the longest predicted peptide
    if args.longest: 
        ep_pred_df = ep_pred_df.sort_values('length').drop_duplicates(['chr', 'pos','gene'], keep='last')
    verboseprint("Kept longest peptide per gene", ep_pred_df.head())
    
    # write modifed df as <qbic_barcode>_modified_ts
    ## replace  NAN values with empty strings
    # ep_pred_df = ep_pred_df.fillna("")

    # drop rows where columns 'chr' or 'pos' are blank
    ep_pred_df = ep_pred_df.dropna(subset = ['chr', 'pos'])  
    
    # sort by columns "QBIC_barcode", then "chr", "length"  
    ep_pred_df  = ep_pred_df.sort_values(["QBIC_barcode", "chr", "length"])
    verboseprint("dimension of ep_pred_df to be written ", ep_pred_df.shape)
    
    outfile_ep_pred_df = qbic_barcode +  "_prediction_result_" + args.suffix + ".tsv"
    fullpath_outfile_ep_pred_df  = os.path.join(args.path_abs_to_seq, dirs, outfile_ep_pred_df)
    ep_pred_df.to_csv(fullpath_outfile_ep_pred_df, sep="\t", index=False)
     
    # compile WB and SB epitopes
    ## get columns of "HLAxx:xx rank"
    ## retrieve regex of "HLAxx:xx", then add HLA count in one of the 4 dicts
    ep_binding_df = ep_pred_df.filter(regex="HLA-[ABC]\*\d{1,3}:\d{1,4} rank")
    # convert ep_binding_df to list of dictionaries (HLAs)
    hlaCols_dicts_dict = ep_binding_df.to_dict(orient="list", into=OrderedDict)
    verboseprint("hlaCols_dicts_dict ", hlaCols_dicts_dict)
    
    # Create 2 lists, SB and WB with of length(dict)
    WB = [0] * len(ep_binding_df.index)
    SB = [0] * len(ep_binding_df.index)

    # set condition of patient
    if sample_status == "metastasis":
        condition = "metastasis"
    else:
        condition = "primaryTumour"
    # verboseprint("sample_status before SB/WB binding ", sample_status)
    # Iterate over each dictionary (HLAs), add +1 based on position
    for hla_keys, hla_ranks in hlaCols_dicts_dict.items():
        hla  = hla_keys.split(" ")[0]
        verboseprint("hla ", hla_keys)
        for count, eachRank in enumerate(hla_ranks):
            if eachRank <= 0.5:
                SB[count] += 1
                sb_key = condition + "_SB_" + hla
                update_dict_counter(ep_binding_hlas_dict, sb_key)
            elif eachRank > 0.5 and eachRank <= 2.0:
                WB[count] += 1
                wb_key = condition + "_WB_" + hla
                update_dict_counter(ep_binding_hlas_dict, wb_key)
        verboseprint("SB list " , SB)
        verboseprint("WB list " , WB)
        verboseprint("ep_bind_hlas_dict ", list(ep_binding_hlas_dict.items()))
    
    # write out SB_list and WB_list into an output table:
    ## columns: QBIC_barcode | patient | condition | sequence | gene | SB | WB
    binding_df = df2 = ep_pred_df.loc[:, ["QBIC_barcode","patient","condition","sequence", "gene", ]]
    binding_df["SB"] = SB # add columns
    binding_df["WB"] = WB  

    outfile_ep_binding_df = qbic_barcode +  "_binding_rank_" + args.suffix + ".tsv"
    fullpath_outfile_ep_binding_df  = os.path.join(args.path_abs_to_seq, dirs, outfile_ep_binding_df)
    binding_df.to_csv(fullpath_outfile_ep_binding_df, sep="\t", index=False)
 
    # create dynamic updating tables (for each sample) for merge_tables and HLA_distribution
    if args.mergeTables:
        # A. Compiling information for common merged tables if args.mergeTables
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
            if cols_ep_pred_df_set != overlap_colNames_set: # more columns in ep_pred_df
                complement_ep_pred_df_set = cols_ep_pred_df_set.difference(overlap_colNames_set)
                verboseprint(" colNames ep_pred_df before dropping columns ", ep_pred_df.columns.tolist()[-7:])
                mod_ep_pred_df = ep_pred_df[ep_pred_df.columns.difference(list(complement_ep_pred_df_set), sort = False)]
                verboseprint(" colNames mod_ep_pred_df after dropping columns ", mod_ep_pred_df.columns.tolist()[-7:])
            if cols_merged_ep_pred_df_set != overlap_colNames_set: # more columns in merged_ep_pred_df
                complement_merged_ep_pred_df_set = cols_merged_ep_pred_df_set.difference(overlap_colNames_set)
                verboseprint(" colNames merged_ep_pred_df before dropping columns ", merged_ep_pred_df.columns.tolist()[-7:])
                mod_merged_ep_pred_df = merged_ep_pred_df[merged_ep_pred_df.columns.difference(list(complement_merged_ep_pred_df_set), sort = False)]
                verboseprint(" colNames mod_merged_ep_pred_df after dropping columns ", mod_merged_ep_pred_df.columns.tolist()[-7:])
            else:
                # the columns in overlap_colNames_set are the same with merged_ep_pred_df and ep_pred_df
                mod_merged_ep_pred_df = merged_ep_pred_df.copy()
                mod_ep_pred_df = ep_pred_df
            
            # append mod_ep_pred_df to mod_merged_ep_pred_df
            # just to be safe, make sure the order of columns in ep_pred_df = merged_ep_df           
            mod_merged_ep_pred_df = mod_merged_ep_pred_df[mod_merged_ep_pred_df.columns.tolist()]
            # mod_ep_pred_df = mod_ep_pred_df[mod_merged_ep_pred_df.columns.tolist()]
            merged_ep_pred_df = pd.concat([mod_merged_ep_pred_df, mod_ep_pred_df])
            verboseprint('Dimenions of concatenated merged_ep_pred_df ', merged_ep_pred_df.shape)
            # reset colNames_list 
            colNames_list = merged_ep_pred_df.columns.tolist()
        
        # B. Compiling information for HLA_distribution
        # write out table HLA distribuition in the form of :
        # Allele | SB | HB
        ## create dicts for SB and WB, respectively
        hla_dist_bind_primaryTumor_dict = {}
        hla_dist_bind_metastasis_dict = {}

        for hla_bind_info, bind_counts in ep_binding_hlas_dict.items():
            verboseprint("hla_bind_info ", hla_bind_info)
            # split hla_bind_info
            hla_condition, bind_strength, hla_coord = hla_bind_info.split("_")
            if hla_condition == "metastasis":
                update_global_hla_counter(hla_dist_bind_metastasis_dict, hla_coord, bind_strength, bind_counts)
            elif hla_condition == "primaryTumour":
                update_global_hla_counter(hla_dist_bind_primaryTumor_dict, hla_coord, bind_strength, bind_counts)
 


# write out merged_table and HLA_distribution_table across samples
if args.mergeTables:
    # A. write out mergedTables
    merged_ep_pred_df = merged_ep_pred_df.fillna("")
    # sort by column "QBIC_barcode", then "chr", "length"
    merged_ep_pred_df  = merged_ep_pred_df.sort_values(["QBIC_barcode", "chr", "length"])
    verboseprint("dimension of merged_ep_pred_df (concat) to be written \n ", merged_ep_pred_df.shape)
    fullpath_outfile_merged_ep_pred_df = os.path.join(args.path_abs_to_seq, "merged_ep_pred_df.tsv")
    merged_ep_pred_df.to_csv(fullpath_outfile_merged_ep_pred_df, sep="\t", index=False)

    # B. write out HLA_distribution table
    verboseprint("hla_dist_bind_metastasis_dict ", hla_dist_bind_metastasis_dict.items())
    verboseprint("hla_dist_bind_primaryTumor_dict ", hla_dist_bind_primaryTumor_dict.items())

    # create DFs from dicts, sort by HLA before writing out
    verboseprint("hla_dist_bind_metastasis_dict, ",  list(hla_dist_bind_metastasis_dict.items()))
    hla_dist_metastasis_df = pd.DataFrame.from_dict(hla_dist_bind_metastasis_dict, orient="index",
                                                   columns = ["SB", "WB"])
    # convert rows (index) to column "HLA", sort HLA before writing
    hla_dist_metastasis_df.reset_index(inplace=True)
    hla_dist_metastasis_df = hla_dist_metastasis_df.rename(columns = {'index' : "HLA"})
    hla_dist_metastasis_df = hla_dist_metastasis_df .sort_values(["HLA"])
    verboseprint('hla_dist_metastasis_df' ,  hla_dist_metastasis_df.head())

    fullpath_outfile_hla_dist_metastasis_df  = os.path.join(args.path_abs_to_seq, "metastasis_SB_WB.tsv")
    hla_dist_metastasis_df.to_csv(fullpath_outfile_hla_dist_metastasis_df, sep="\t", index=False)

    # repeat for primaryTumour
    hla_dist_primaryTumour_df = pd.DataFrame.from_dict(hla_dist_bind_primaryTumor_dict, orient="index",
                                                 columns = ["SB", "WB"])
    hla_dist_primaryTumour_df.reset_index(inplace=True)
    hla_dist_primaryTumour_df = hla_dist_primaryTumour_df.rename(columns = {'index' : "HLA"})
    hla_dist_primaryTumour_df = hla_dist_primaryTumour_df .sort_values(["HLA"])
    verboseprint('hla_dist_primaryTumour_df' ,  hla_dist_primaryTumour_df.head())

    fullpath_outfile_hla_dist_primaryTumour_df = os.path.join(args.path_abs_to_seq, "primaryTumour_SB_WB.tsv")
    hla_dist_primaryTumour_df.to_csv(fullpath_outfile_hla_dist_primaryTumour_df, sep="\t", index=False)









    



                                                                       
                                                                       

    





