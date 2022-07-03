#!/usr/bin/env python
# Jun-Hoe, Lee (2020)
# spin-off code snippet to check number of min sequences before copying to another folder
# requires biopython package
# read a fasta file, and based on the given minimum seq, write out file in single line
# usage: 
## python filterMinSpecies.py [inputSeqFile] [outputSeqFile] [totalSpecies] -minNumSpecPer -verbose

import os
import sys
import subprocess
import argparse
import configparser

from Bio import SeqIO

#####################################################
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


parser = errorDisplayParser(description='filterMinSpecies: Read a fasta file, check number of seqs, \
                                         then output if passed min number of sequences ')
parser.add_argument('inputSeqFile', action='store' , help='input fasta file')
parser.add_argument('outputSeqFile', action='store',
                     help='output fasta file (single line) if passed min sequences')
parser.add_argument('totalSpecies', action="store", type=int, help='total number of species in all gene sets')
parser.add_argument('-minNumSpecPer', action="store", type=int, default=80, help='minimum number of species (percent of totalSpecies, default = 80 percent)')
parser.add_argument('-verbose', action='store_true', default=False, help='turns on verbose mode. Usage: -verbose')  # verbose flag

args = parser.parse_args()

# function to turn on verbose mode
if args.verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print("-v ", arg, ":", args)
else:
    verboseprint = lambda *a: None  # do-nothing function

#########################################################
# A. read input fasta files


seqData = list(SeqIO.parse(args.inputSeqFile, "fasta"))
verboseprint("read input, seq 1 ", seqData[0].id, seqData[0].seq)
verboseprint("len seqData (number of seqs) ", len(seqData))

## check number of min seqs
minSpecies = round(args.minNumSpecPer * args.totalSpecies / 100)
verboseprint("minSpecies ", minSpecies)

# write output if exceed filterMinSpecies

newLine = ''
if len(seqData) >= minSpecies:
    with open(args.outputSeqFile, 'w') as writeOutputFile:
        for eachSeq in seqData:
            writeLine = ">" + str(eachSeq.id) + "\n" + str(eachSeq.seq)
            writeOutLine = newLine + writeLine
            newLine = '\n'
            writeOutputFile.write(writeOutLine)
    writeOutputFile.close()
else:
    print(args.inputSeqFile + "did not exceed min number of seqs ", len(seqData))