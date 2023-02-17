#!/usr/bin/env python
# Jun-Hoe, Lee (2022)
# script to create sample sheet for nf-corek pipeline, specifically project QHPRG
# modified from script_sarek_samplesheet.py - now to reannotate VCF files from IMGAG
# input: i. list (tsv) of QBIC_barcode, sample_grouping (Fxxxxx), and sample_type
# ii. path to directories (QBIC_barcodes) containing the VCF files
# iii. output sample sheet name 
# output: sample sheet file (csv) with 3 columns: patient,sample,vcf(file location)
# usage: 
# python sarek_create_samplesheet.py [sample_info_list] [full path to samples] [out_samplesheet] -verbose


import io
import os
import pathlib
import sys
import subprocess
import argparse
import configparser
import pandas as pd
import re
import mmap
import gzip


from io import StringIO
from collections import OrderedDict



######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Create samplesheet for nf-core sarek pipeline, specifically project QHPRG')
parser.add_argument('input_sample_metadata', action="store", type=argparse.FileType('r'),
                     help="tsv: QBIC_barcodes\'\t\'sample_grouping\'\t\'sample_type")
# parser.add_argument('tumour_type', action="store",
#                      help="primary tumor 'pt' or metastasis 'ms' ")                     
parser.add_argument('path_abs_to_seq', action='store',
                     help="absolute path to main directory that contains dirs of VCF files")
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
# function to create a copy of VCF files, replacing 'CSQ' with 'CSQ_IMGAG' 
def modify_CSQ_in_VCF(old_filename, new_filename, old_string, new_string):
    with gzip.open(old_filename,"rt") as rf:
        newText=rf.read().replace(old_string, new_string)

    with gzip.open(new_filename, "wt") as wf:
        wf.write(newText)


# function to name and rename old/new VCF files
def generate_old_new_VCF_filename(dir_path, qbic_sample_barcode, vc_sample): 
    old_vc_file_name = vc_sample + "_var_annotated.vcf.gz"
    path_old_VCF_file = os.path.join(dir_path, qbic_sample_barcode, old_vc_file_name)
    verboseprint("full_path_file VCF ", path_old_VCF_file)

    new_VCF_file = vc_sample + "_var_annotated_updated_To_CSQ_IMGAG.vcf.gz"
    path_new_VCF_file = os.path.join(dir_path, qbic_sample_barcode, new_VCF_file)

    return path_old_VCF_file, path_new_VCF_file


# function to check whether VCF file headers match
def check_VCF_filenames_in_headers(input_VCF_file, VCF_filename):
    with gzip.open(input_VCF_file, "rb", 0) as f:
        datafile = f.read()             
        if str.encode(VCF_filename) in datafile:
            verboseprint("found VCF_filename ", VCF_filename)
            return True
        else:
            verboseprint("couldn't find VCF_filename ", VCF_filename)
            return False
           

def write_full_line(mod_VCF_dict):
    # from patient_modVCF_dict, retrieve all necessary info to generate a tabbed line.
    # e.g. below
    # patient   sample  vcf
    # 3 FO13399x01_01    xx
    output_line_list = []
     
    for patient_num, patient_vcfs in mod_VCF_dict.items():
        for each_vcf in patient_vcfs:
        # unpack the blood, tumour, (and) metastasis samples
            VCF_filename = os.path.basename(each_vcf)
            VCF_sample_name = VCF_filename.removesuffix("_var_annotated_updated_To_CSQ_IMGAG.vcf.gz")
            output_single_list = [patient_num, VCF_sample_name, each_vcf]
            single_line = "\t".join(output_single_list) 
            output_line_list.append(single_line)
        
    verboseprint("output_line_list", output_line_list[:1])
    
    return output_line_list


#########################################################
# A. read in input_sample_metadata, store as dictionary

## example of input_sample_metadata,
# QBiC Code	Organism ID	Analyte ID	Condition: sample_type
# QHPRG586AJ	3	FO13399x01_01	blood
# QHPRG587AR	3	FO13399x02_01	primary tumor
# QHPRG588A1	7	FO13518x03_01	blood
# QHPRG589A9	7	FO13518x02_01	primary tumor


readInputFileNames = args.input_sample_metadata.readlines()[1:] # skip header
sample_metadata_dict = OrderedDict()
for line in readInputFileNames:
    line = line.strip() # remove \n
    splitLine = line.split("\t")
    sample_ID, sample_patient, sample_patient_analyte, sample_condition = splitLine
    if sample_condition == "blood":
        sample_status_label = "bl"
    else:
        sample_status_label = "tu"
    
    sample_patient_status = sample_status_label + "," + sample_ID + "," + sample_patient_analyte
    if sample_patient not in sample_metadata_dict:
        sample_metadata_dict[sample_patient] = [sample_patient_status]
    else:
        sample_metadata_dict[sample_patient].append(sample_patient_status)

verboseprint("sample_metadata_dict ", list(sample_metadata_dict.items())[:3])
## e.g. {3: ["bl, QHPRG586AJ, FO13399x01_01" , "pt/ms, QHPRG586AJ, FO13399x02_01"]...}


# B. walk through the directory of VCF files from IMGAG, modify 'CSQ' field
try:
    dir_list = os.listdir(args.path_abs_to_seq)
    verboseprint("dir_list ", dir_list)
    
except:
    print("Incorrect directory path provided. ")
    sys.exit()

# check-match my own generated file_names with existing files, before modify vcf file
# each samples are named with the convention, e.g.:
## QHPRG586AJ (blood)
## QHPRG587AR_QHPRG586AJ (primary tumour vs blood)
## Qxxxx_QHPRG586AJ (metastasis vs blood)

# new dict to store patient_num: modified_vcf_files
patient_modVCF_dict = OrderedDict()

for patient_num, patient_metadata in sample_metadata_dict.items():
    # sort the nested lists by status, i.e. "bl" should be first
    sorted_sample_list = sorted(patient_metadata, key=lambda x: x[0])
    verboseprint("sorted_sample_list ", sorted_sample_list)
    # deal with with blood sample first, label each item there
    status_blood, qbic_barcode_blood, vc_sample_blood = sorted_sample_list[0].split(",")
    old_vcf_file_blood, new_vcf_file_blood = generate_old_new_VCF_filename(args.path_abs_to_seq, qbic_barcode_blood, vc_sample_blood)
    verboseprint("old_vcf_file_blood ", old_vcf_file_blood)
    verboseprint("new_vcf_file_blood ", new_vcf_file_blood)

    if os.path.exists(old_vcf_file_blood) and check_VCF_filenames_in_headers(old_vcf_file_blood, vc_sample_blood):
        # check that the VCF file names, Fxxxx-Fxxx is in the VCF file:
        modify_CSQ_in_VCF(old_vcf_file_blood, new_vcf_file_blood, "CSQ", "CSQ_IMGAG")
        patient_modVCF_dict[patient_num] = [new_vcf_file_blood]
    
    # now modify the primary tumour, and (if exist) metastasis samples
    for each_sample in sorted_sample_list[1:]: # skip blood sample
        status, qbic_barcode, vc_sample = each_sample.split(",")
        qbic_barcode_tumour = qbic_barcode + "_" + qbic_barcode_blood
        vc_sample_tumour = vc_sample + "-" + vc_sample_blood
        old_vcf_file, new_vcf_file = generate_old_new_VCF_filename(args.path_abs_to_seq, qbic_barcode_tumour, vc_sample_tumour)
        if os.path.exists(old_vcf_file) and check_VCF_filenames_in_headers(old_vcf_file, vc_sample_tumour) :
            modify_CSQ_in_VCF(old_vcf_file, new_vcf_file, "CSQ", "CSQ_IMGAG")
            patient_modVCF_dict[patient_num].append(new_vcf_file)

verboseprint("patient_modVCF_dict ", list(patient_modVCF_dict.items())[:3])


# C. create samplesheet, in tabbed lines

## example of output_sample_sheet
# patient   sample  vcf
# 3 FO13399x01_01    xx
# 3 FO13399x02_01   yy 

write_lines_list = write_full_line(patient_modVCF_dict)


# D. import write_lines_list into pandas,export as csv
header_names = ["patient","sample","vcf"]

# samplesheet_df = pd.read_csv(StringIO(write_lines_list), sep='\t', header=header_names)
samplesheet_df = pd.read_csv(io.StringIO('\n'.join(write_lines_list)), names = header_names, sep="\t")

verboseprint("samplesheet_df head ", samplesheet_df.head())

# write to csv
samplesheet_df.to_csv(args.output_samplesheet, index=False)




