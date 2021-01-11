#/usr/bin/python3
# Jun-Hoe, Lee (2019)
# convert multiline fasta to a single file
# usage: python format_fasta_oneline.py [input file] [- outputfile]

import os
import sys
import subprocess
import argparse
import configparser



######################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = errorDisplayParser(description='Converts a multline fasta file into a single line fasta')
parser.add_argument('inputSeqFile', action='store',
                     help='input file containing a multiline fasta')
parser.add_argument('-output', type=argparse.FileType('w'),
                     help='output file converted to a single line fasta')
parser.add_argument('-verbose', action='store_true', default=False, help='turns on verbose mode. Usage: -verbose')  # verbose flag


args = parser.parse_args()


# function to turn on verbose mode
if args.verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print("-v ", arg, ":", args)
else:
    verboseprint = lambda *a: None  # do-nothing function

#########################################################
# run command

try:
    callConvertFormat = ("awk \'{if(NR==1) {print $0} else {if($0 ~ /^>/) {print \"\\n\"$0} else {printf $0}}}\' " + args.inputSeqFile)
    verboseprint("callConvertFormat", callConvertFormat)
    if args.output:
        procConvertFormat = subprocess.Popen(callConvertFormat, stdout= subprocess.PIPE, shell = True) # pipe the output from stdout
        convertStdOut = procConvertFormat.stdout.read() # the stdOut is a fasta file in very long single line
        decodeOutput = convertStdOut.decode('utf-8') # decode from bytes to string
        splitConvertStdOut = decodeOutput.split("\n") # split into list using "\n", which are found between the species name and sequence
    else:
#        procConvertFormat = subprocess.run(callConvertFormat, shell =True)
        subprocess.run(callConvertFormat, shell=True)
except:
    print("failed to convert to single line fasta")

# print out output file

if args.output:
    newLine = ''
    for line in splitConvertStdOut:
        line = line.strip() # remove any possible starting and trailing 
        writeLine = newLine + line
        verboseprint("write line to file ", line)
        args.output.write(writeLine)
        newLine = '\n'
