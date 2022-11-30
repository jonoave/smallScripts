#!/usr/bin/env python
# Jun-Hoe, Lee (2022)
# script to create sample sheet for nf-corek pipeline, specifically project QHPRG
# input: i. list (tsv) of QBIC_barcode, sample_grouping, and sample_type
# ii. path to directories (QBIC_barcodes) containing the fastq.gz
# iv. full path to samples (str)  
# output: sample sheet file (csv)
# usage: 
# python sarek_create_samplesheet.py [sample_info_list] [path to directory] [full path to samples] [out_samplesheet] -verbose


import io
import os
import pathlib
import sys
import subprocess
import argparse
import configparser
import pandas as pd
import re
from io import StringIO


######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Create samplesheet for nf-core sarek pipeline, specifically project QHPRG')
parser.add_argument('input_sample_metadata', action="store", type=argparse.FileType('r'),
                     help="tsv: QBIC_barcodes\'\t\'sample_grouping\'\t\'sample_type")
parser.add_argument('path_abs_to_seq', action='store',
                     help="absolute path to main directory that contains dirs of seqs")
parser.add_argument('output_samplesheet', action='store',
                     help="output samplesheet file name")
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

def tag_r1_r2(qbic_barcode, lane_r1, lane_r2, lane_x, list_of_files):
    # function to extract the corresponding R1,R2 files given a lane number (L00x)
    # returns a list of [L00x_(L00x_R1)_(L00x_R2)]

    lane_R1_file = [r1 for r1 in list_of_files if lane_r1 in r1][0]
    lane_R2_file = [r2 for r2 in list_of_files if lane_r2 in r2][0]
    # append location (QBIC barcode + "/") prefix to lane_R1/R2
    lane_R1_file = qbic_barcode + "/" + lane_R1_file 
    lane_R2_file = qbic_barcode + "/" + lane_R2_file
    list_one_lane_r1_r2 = [lane_x, lane_R1_file, lane_R2_file]
    verboseprint("list_one_lane_r1_r2 ", list_one_lane_r1_r2)
    
    return list_one_lane_r1_r2


def parse_read_lane_files(list_of_files):
    # function to parse list of files based on R1/2 and L00x
    # identify if or how many L00x among samples in list_of_files
    # returns list of list: [[L00x,L00x_R1,L00x_R2], [L00y,L00y_R1,L00y_R2]...]

    find_L00x = re.compile(r"(Q\w{9})_(L00\d)")
    list_L00x = []


    lane_r1_r2_list = []
      
    # check for presence of lanes
    for samples_files in list_of_files:
        match_L00x = find_L00x.match(samples_files)
        if match_L00x:
            list_L00x.append(match_L00x.group(2))

    # if there are different lanes:
    if list_L00x: # not empty list 
        sample_ID_loc = match_L00x.group(1)
        # convert list_L00x to set, then sort by ascending
        set_L00x = sorted(set(list_L00x))
        verboseprint("set_L00x ", set_L00x)
        
        # start creating list of [Lx_(Lx_R1_file_name)_(Lx_R2_file_name),]
        ## pull out the respective R1 and R2 files
        for lanes in set_L00x:
            lane_R1 = lanes + "_R1"
            lane_R2 = lanes + "_R2"
            list_lane_r1_r2 = tag_r1_r2(sample_ID_loc, lane_R1, lane_R2, lanes, list_of_files)
            lane_r1_r2_list.append(list_lane_r1_r2)
    # if there are no lanes (no L00x), then just tag R1, R2 samples in one line
    else:
        sample_ID_loc = list_of_files[0][:10]
        list_lane_r1_r2 = tag_r1_r2(sample_ID_loc, "R1", "R2", "L001", list_of_files)
        lane_r1_r2_list.append(list_lane_r1_r2)
    
    return lane_r1_r2_list


def write_full_line(qbic_id_sample, metadata_dict, file_list):
    # from qbic_ID_sample, retrieve all necessary info to generate a tabbed line
    # in samplesheet in the order of "patient,status,sample,lane,fastq_1,fastq_2"

    sample_info =  metadata_dict[qbic_id_sample]
    patient_id = sample_info.split("_")[0]
    sample_status = sample_info.split("_")[1]
    sample_id = qbic_id_sample

    # parse the list of files and group them by lane, R1/R2
    all_lane_r1_r2_files_list = parse_read_lane_files(file_list)
    
    # list to store each line of samplesheet
    output_line_list = []

    for each_lane in all_lane_r1_r2_files_list:
        lane, r1, r2 = each_lane
        output_single_line = [patient_id, sample_status, sample_id, lane, r1, r2]
        single_line = "\t".join(output_single_line) 
        output_line_list.append(single_line)

    verboseprint("output_line_list", output_line_list[:1])

    return output_line_list



#########################################################
# A. read in input_sample_metadata, store as dictionary

## example of input_sample_metadata,
# QHPRG586AJ  3   blood
# QHPRG587AR  3   primary tumor
# QHPRG594AA  10  metastasis
# QHPRG595AI  10  primary tumor
# QHPRG596AQ  10  blood

# sample status: {"blood":"0", "primary tumor":"1", "metastasis":"1"} 

readInputFileNames = args.input_sample_metadata.readlines()[1:] # skip header
sample_metadata_dict = {}
for line in readInputFileNames:
    line = line.strip() # remove \n
    splitLine = line.split("\t")
    sample_ID, sample_patient, sample_status = splitLine
    if sample_status == "blood":
        sample_status_label = "0"
    else:
        sample_status_label = "1"
    sample_patient_status = sample_patient + "_" + sample_status_label
    sample_metadata_dict[sample_ID] = sample_patient_status

verboseprint("sample_metadata_dict ", list(sample_metadata_dict.items())[:3])
## e.g. {"QHPRG586AJ": "3_0", "QHPRG586AJ": "3_1"]...}

# B. walk through the directory of NGS to look for fastq files
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

for dirs in dir_list:
    if dirs in sample_metadata_dict:
        fullpath = os.path.join(args.path_abs_to_seq, dirs)
        # only append files ending with fastq.qz, and excluding those with fastq.gz.crc32
        # get list of files 
        files_list = os.listdir(fullpath)
        
        for files_name in files_list:
            if "fastq.gz" in files_name and "crc32" not in files_name:
                if dirs not in sample_seqFiles_dict:
                    sample_seqFiles_dict[dirs] = [files_name]
                else:
                    sample_seqFiles_dict[dirs].append(files_name)


verboseprint("sample_seqFiles_dict  ", list(sample_seqFiles_dict .items())[:3])

# remove entries (samples) that are lacking either R1/R2
verboseprint("number of entries (samples) in sample_seqFiles_dict: ", len(sample_seqFiles_dict))
sample_seqFiles_paired_dict = {k:v for k, v in sample_seqFiles_dict.items() if len(v) % 2 == 0 }
verboseprint("number of entries (samples) in sample_seqFiles_paired : ", len(sample_seqFiles_paired_dict))


# C. create sample sheet
## sample sheet is a csv file with the following columns:
## patient,status,sample,lane,fastq_1,fastq_2

## first create write_lines_list (tabbed) that stores each line of the sample sheet
## second, import write_lines into a pandas dataframe. Then sort by patient, status, R1, lane etc
## third, append complete R1,R2 pathway string to R1,R2 column
## finally, export as csv file
           
## create write_lines_list
write_lines_list = []
  
for qbic_id, qbic_files in sample_seqFiles_paired_dict.items():
    lines_list = write_full_line(qbic_id, sample_metadata_dict, qbic_files)
    # this returns a list of lists, could be only 1 element if there is only a single lane 
    write_lines_list.extend(lines_list)

verboseprint("write_lines_list ", write_lines_list[:5])

# D. import write_lines_list into pandas to sort by columns
header_names = ["patient","status","sample","lane","fastq_1","fastq_2"]

# samplesheet_df = pd.read_csv(StringIO(write_lines_list), sep='\t', header=header_names)
samplesheet_df = pd.read_csv(io.StringIO('\n'.join(write_lines_list)), names = header_names, sep="\t")

verboseprint("samplesheet_df head ", samplesheet_df.head())

# sort by columns "patient", "status" then by "lane"
samplesheet_sort_df = samplesheet_df.sort_values(["patient", "status", "lane"])

verboseprint("samplesheet_sort_df head ", samplesheet_sort_df.head())

# add absolute path to columns "fastq1" and "fastq2"
samplesheet_sort_df["fastq_1"] = args.path_abs_to_seq + "/" + samplesheet_sort_df["fastq_1"].astype(str)
samplesheet_sort_df["fastq_2"] = args.path_abs_to_seq + "/" + samplesheet_sort_df["fastq_2"].astype(str)

verboseprint("samplesheet_sort_df head, added data path ", samplesheet_sort_df.head())

# write to csv
samplesheet_sort_df.to_csv(args.output_samplesheet, index=False)




