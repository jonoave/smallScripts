#!/usr/bin/env python
# Jun-Hoe, Lee (2024)
# script to find files recursively deep in folders and move them to <output> dir
# input: i. input folder to traverse recursively
# ii. output folder to copy files to 
# usage:
# python findFilesRecursivelyAndMove.py [input_folder]] [file_prefix] [outputFolder'] -verbose
# WARNING! Make backup before executing script to rename!

import os
import shutil
import sys
import argparse


######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Traverse directories recursively, from input_folder and copy all files with \
                            file_prefix to output_folder ')
parser.add_argument('input_folder', action='store',
                     help="absolute path to main starting directory")
parser.add_argument('file_prefix', action='store',
                     help='constant string that appears in all files')
parser.add_argument('output_folder', action='store',
                     help='path to directory')
parser.add_argument('-verbose', action='store_true', default=False,
                    help='turns on verbose mode. Usage: -verbose')  # verbose flag


args = parser.parse_args()

# function to turn on verbose mode
if args.verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print("-v ", arg, ":", args)
else:
    verboseprint = lambda *a: None  # do-nothing function

# function to traverse files 
def walkfs(startdir, findfile):
    list_of_files = []
    verboseprint("created list of files " , list_of_files)
    for root, dirs, files in os.walk(startdir):
        verboseprint("root ", root)
        verboseprint("dirs ", dirs)
        verboseprint("files", files)
        if files: # files list is not empty
            file_list = [x for x in files if findfile in x]
            # add fullpath to each file name
            file_list_root = [os.path.join(root, x) for x in file_list]
        list_of_files += file_list_root

    return list_of_files


def createDestinationFolder(destination):
    isExist = os.path.exists(destination)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(destination)
        verboseprint("The new directory is created! \n", destination)


#########################################################

## using template from here: 
## https://stackoverflow.com/questions/39879638/recursively-go-through-all-directories-until-you-find-a-certain-file-in-python


# get file names recursively
files_list = walkfs(args.input_folder, args.file_prefix)
verboseprint("files_list ", files_list)

# now add root directory to each item in list
root_path = os.getcwd()
file_list_full_root = [os.path.join(root_path, x) for x in files_list]
verboseprint("files_list " , file_list_full_root[:2])


# now copy files 

output_path = os.path.join(root_path,args.output_folder)
verboseprint(output_path)
createDestinationFolder(output_path )

for each_file in file_list_full_root:
     shutil.copy(each_file, output_path)

