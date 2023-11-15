#!/usr/bin/env python
# Jun-Hoe, Lee (2023)
# Designed for new batch of data transfered from IMGAG (October 2023)
# script to parse through list of files, compare them to samplesheet":
# i. rename Fxxxx files to QBIC barcodes
# ii. remove 01_Sxxx from sample names, 
# iii. change 001.1 type labels to L001_R1 
# iv. group all reads of one sample in a sample folder 
# usage:
# python renameLabelMoveFiles.py [samplesheet] [source directory] [target directory] -verbose
# WARNING! Make backup before executing script to rename!


import os
import pathlib
import sys
import subprocess
import argparse
import configparser
import re

from shutil import copyfile


######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Rename files and group reads by samples in target directory')
parser.add_argument('input_sample_metadata', action="store", type=argparse.FileType('r'),
                     help="tsv: QBIC_barcodes\'\t\'sample_grouping\'\t\'sample_type")
parser.add_argument('path_abs_to_seq', action='store',
                     help="absolute path to main directory that contains all files")
parser.add_argument('output_path', action='store',
                     help="absolute path to create individual sample folders")
parser.add_argument('-verbose', action='store_true', default=False,
                    help="turns on verbose mode. Usage: -verbose")  # verbose flag



args = parser.parse_args()
# function to turn on verbose mode
if args.verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print("-v ", arg, ":", args)
else:
    verboseprint = lambda *a: None  # do-nothing function


########################################################
# function to check if regex is present
def is_match(regex, text):
    pattern = re.compile(regex)
    match_pattern = pattern.match(text)
    return match_pattern

def createDestinationFolder(destination):
    isExist = os.path.exists(destination)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(destination)
        verboseprint("The new directory is created! \n", destination)


#########################################################
# A. read in input_sample_metadata, store as dictionary

## example of input_sample_metadata,
# QBiC Code	Organism ID	Analyte ID	Condition: sample_type
# QHPRG586AJ	3	FO13399x01_01	blood
# QHPRG587AR	3	FO13399x02_01	primary tumor
# QHPRG588A1	7	FO13518x03_01	blood
# QHPRG589A9	7	FO13518x02_01	primary tumor


# create 2 dictionaries
# analyte_to_qbic_dict: e.g FO13399x01_01:QHPRG586AJ
# qbic_to_metadata_dict: e.g. QHPRG586AJ: ["3", "bl", "FO13399x01_01"]

readInputFileNames = args.input_sample_metadata.readlines()[1:] # skip header

analyte_to_qbic_dict = {}
qbic_to_metadata_dict = {}

for line in readInputFileNames:
    line = line.strip() # remove \n
    splitLine = line.split("\t")
    sample_ID, sample_patient, sample_patient_analyte, sample_condition = splitLine
    if sample_condition == "blood":
        sample_status_label = "bl"
    else:
        sample_status_label = "tu"
    
    # create entry in analyte_to_qbic_dict:
    analyte_to_qbic_dict[sample_patient_analyte] = sample_ID
    
    # create entry in qbic_to_metadata_dict
    sample_patient_metadata = [sample_patient, sample_status_label, sample_patient_analyte]
    qbic_to_metadata_dict[sample_ID] = sample_patient_metadata


verboseprint("analyte_to_qbic_dict ", list(analyte_to_qbic_dict.items())[:4])
verboseprint("qbic_to_metadata_dict ", list(qbic_to_metadata_dict.items())[:4])



# B. walk through files, modify and move file
# get list of files
try:
    files_list = os.listdir(args.path_abs_to_seq)
    verboseprint("files_list ", files_list)
    
except:
    print("Incorrect directory path provided. ")
    sys.exit() 

re_analyte = "FO\d+"

# dict to convert number lane/reads to letters
lane_reads_dict = {"001.1": "L001_R1_001",
                   "001.2": "L001_R2_001",
                   "002.1": "L002_R1_001",
                   "002.2": "L002_R2_001"}
verboseprint("lane_reads_dict ", list(lane_reads_dict.items()))

for fileNames in files_list:
    # Check for analyte ID, or else use QBIC barcode
    check_match = is_match(re_analyte, fileNames)
    verboseprint("check_match ", check_match)
    if check_match: # e.g. Fxxxxxx_01_S43_L001_R2_001.fastq.gz
        # get corresponding QBIC barcode to analyte_ID
        fileNames_split= fileNames.split("_")
        analyte_ID = fileNames_split[0] + "_" + fileNames_split[1]
        verboseprint("analyte_ID of check_match ", analyte_ID)
        qbic_code = analyte_to_qbic_dict.get(analyte_ID)
        verboseprint("qbic_code if check_match found: ", qbic_code)
        #  now check for lane and read
        # change file name to QHPRGxxxx_L00x_Rx_001.fastq.gz
        # newFileName = fileNames.replace(check_match.group(0), qbic_code)
        newFileName = qbic_code + "_" + '_'.join(fileNames_split[-3:]) 
        verboseprint("newFileName if check_match ", newFileName)
    else: # e.g. QHPRGxxxx_normal_001.1.fastq.gz
        # split file name look for "tumor" or "normal"
        splitFilename = fileNames.split("_")
        qbic_code, sample_status, read_lane = splitFilename
        read_lane = read_lane.replace(".fastq.gz", "")
        verboseprint("qbic_code, sample_status, read_lane from splitFileName: \n ", \
                     qbic_code, sample_status, read_lane)
        verboseprint("qbic_code if check_match not found: ", qbic_code)
        get_status = qbic_to_metadata_dict[qbic_code][1] # 2nd item in list
        verboseprint("get_status of qbic_code ", get_status)
        if (get_status == "bl" and sample_status == "normal") \
        or (get_status == "tu" and sample_status == "tumor"):
            newLane = lane_reads_dict.get(read_lane)
            newFileName = qbic_code + "_" + newLane + ".fastq.gz"
            verboseprint("newFileName of QBIC: ", newFileName)
        else:
            continue 

    newFileDestination = os.path.join(args.output_path, qbic_code)
    createDestinationFolder(newFileDestination)
    verboseprint("newFileDestination \n", newFileDestination)
    verboseprint("newFileName ", newFileName)
    newFileDestinationCopy = os.path.join(newFileDestination, newFileName)
    # copy file to target location
    copyfile(os.path.join(args.path_abs_to_seq,fileNames), newFileDestinationCopy)  


