#!/usr/bin/env python
# Jun-Hoe, Lee (2022)
# script to create sample sheet for nf-corek pipeline/epitopeprediction, specifically project QHPRG
# also create a copy of alleles from hlatyping output, and modify it so it can be parsed by the ep pipeline
# input: i. list (tsv) of QBIC_barcode, sample_grouping (Fxxxxx), and sample_type
# ii. path to directory of hlatyping output
# iii. path to directory of (re-)annotated VCF files from sarek output
# output: sample sheet file (csv) with 4 columns: sample,alleles,mhc_class,filename
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

parser = errorDisplayParser(description='Create samplesheet for nf-core epitopeprediction pipeline, specifically project QHPRG')
parser.add_argument('input_sample_metadata', action="store", type=argparse.FileType('r'),
                     help="tsv: QBIC_barcodes\'\t\'sample_grouping\'\t\'sample_type")
parser.add_argument('path_abs_to_alleles', action='store',
                     help="absolute path to main directory that contains dirs of hlatyping pipeline")
parser.add_argument('path_abs_to_annotated_vcf', action='store',
                     help="absolute path to main directory that contains dirs of annotated VCFS from sarek pipeline")
parser.add_argument('output_samplesheet', action='store',
                     help="output samplesheet file name")
parser.add_argument('-vcf_mod_suffix', action='store_true', default="_var_annotated.vcf.gz",
                     help:"suffix for VCF files, if reannotated or named differently from default output of VEP files from sarek")
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

# function to get name of old allele file and generate new names.
def generate_old_new_allele_filenames(dir_path_alleles, qbic_barcode)
    fullpath_alleles = os.path.join(dir_path_alleles, qbic_barcode )
    old_allele_file_name = qbic_barcode + "_result.tsv"
    path_old_allele_file_name  = os.path.join(fullpath_alleles, old_allele_file_name)

    new_allele_file_name = qbic_barcode + "_result.txt"
    path_new_allele_file_name  = os.path.join(fullpath_alleles , new_allele_file_name)

    return path_old_allele_file_name, path_new_allele_file_name


# function to create a copy of alleles from hlatyping output,
## and modify it into single line of alleles, split by ";"
def modify_alleles_in_hlatyping(old_filename, new_filename):
   # read old_allele_file_name, make changes
    with open(old_filename,"rt") as alleles_rf:
        read_old_alleles = alleles_rf.readlines()[1] # skip first line, only read second line
        second_line_split = read_old_alleles.split("\t")
        alleles_list = [x for x in second_line_split if "*" in x] # allele coordinates have "*"

    # open new_allele_file_name
    with open(new_filename, "wt") as alleles_wf:
        # write alleles in a single line, separated by colon
        write_alleles_line = ";".join(alleles_list) 
        alleles_wf.write(write_alleles_line)

# function to write out tabbed lines in list for samplesheet
def write_full_line(mod_VCF_dict):
    # from mod_VCF_dict, retrieve all necessary info to generate a tabbed line # e.g. below
    # sample    alleles mhc_class   filename
    # 3 QHPRG586AJ_FO13399x01_01    I   [annotated vcf file from sarek location]
    # 3 QHPRG5867J_FO13399x02_01    I   ..

    output_line_list = []

    for qbicBarcode, vcf_metadata in mod_VCF_dict.items():
        vcf_files, allele_barcode = vcf_metadata
        # generate vcf_file names, check whether exists:
        vcf_file_dir = qbic_barcode + "_" + vcf_files
        vcf_file_name = vcf_files + args.vcf_mod_suffix
        full_vcf_file_path = os.path.join(args.path_abs_to_annotated_vcf, vcf_file_dir,vcf_file_name)
        verboseprint("full_vcf_file_path (annotated) ", full_vcf_file_path)
        if not os.path.exists(full_vcf_file_path):
            print("%s not found in annotated vcf files " %full_vcf_file_path)
            break
        # generate allele file name, check whether exists:
        allele_file_name = allele_barcode + "_result.txt"
        full_allele_file_path = os.path.join(args.path_abs_to_alleles, allele_barcode, allele_file_name)
        verboseprint("full_allele_file_path " , full_allele_file_path)
        if not os.path.exist(full_allele_file_path):
            print("%s not found in modified alleles " %full_allele_file_path)
            break
        # add information to line list
        output_single_list = [vcf_file_dir, full_allele_file_path, "I", full_vcf_file_path]
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

## purpose is to generate filenames manually, and check against the filenames in part C

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
## e.g. {3: ["bl, QHPRG586AJ, FO13399x01_01" , "tu, QHPRG586AJ, FO13399x02_01"]...}

## create dict of tumour only samples, that will be executed by epitopeprediction pipeline
tumour_qbic_vcf_dict = OrderedDict()

for patient, patient_metadata in sample_metadata_dict.items():
    ## sort the nested lists by status, i.e. "bl" should be first
    sorted_sample_list = sorted(patient_metadata, key=lambda x: x[0])
    qbic_barcode_blood, status_blood, vc_sample_blood = sorted_sample_list[0].split(",")
    ## now get info from "tu" samples
    for each_sample in sorted_sample_list[1:]: # skip blood sample
        status_tumour, qbic_barcode, vc_sample = each_sample.split(",")
        vc_sample_tumour = vc_sample + "-" + vc_sample_blood
        # add qbic_barcode_blood, which carries the hlatyping output of the germline
        tumour_qbic_vcf_dict[qbic_barcode] = [vc_sample_tumour, qbic_barcode_blood]
        

verboseprint("tumour_qbic_vcf_dict ", list(tumour_qbic_vcf_dict.items())[:3])
## e.g. {QHPRG587AJ: [FO13399x02_01_FO13399x01_01, QHPRG586AI ], ...}


# B. walk through the directory of alleles (output from hlatyping)
## create a copy and modify so all alleles are listed in one line, separated by colon
try:
    alleles_dir_list = os.listdir(args.path_abs_to_alleles)
    verboseprint("dir_list ", alleles_dir_list)
except:
    print("Incorrect path to directory of alleles (output of hlatyping pipeline) provided. ")
    sys.exit()


## copy allele files and modify 
for allele_dir in alleles_dir_list: # allele_dir is QBIC_barcode
    old_allele_fileName, new_allele_file_name = generate_old_new_allele_filenames(args.path_abs_to_alleles, allele_dir )
    modify_alleles_in_hlatyping(old_allele_fileName, new_allele_file_name ) 


# C. compile all information into a list for writing samplesheet

## example of output_sample_sheet
# sample    alleles mhc_class   filename
# 3 QHPRG586AJ_FO13399x01_01    I   [annotated vcf file from sarek location]
# 3 QHPRG5867J_FO13399x02_01    I   ..

try:
    vcf_dir_list = os.listdir(path_abs_to_annotated_vcf)
    verboseprint("dir_list ", vcf_dir_list)
except:
    print("Incorrect path to directory of annotated vcfs (output of sarek pipeline) provided. ")
    sys.exit()

write_lines_list = write_full_line(tumour_qbic_vcf_dict)


# D. import write_lines_list into pandas,export as csv
header_names = ["sample","alleles","mhc_class","filename"]

# samplesheet_df = pd.read_csv(StringIO(write_lines_list), sep='\t', header=header_names)
samplesheet_df = pd.read_csv(io.StringIO('\n'.join(write_lines_list)), names = header_names, sep="\t")

verboseprint("samplesheet_df head ", samplesheet_df.head())

## write to csv
samplesheet_df.to_csv(args.output_samplesheet, index=False)




