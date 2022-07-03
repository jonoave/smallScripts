#!/usr/bin/python
# Jun Hoe, Lee (2013)
# parse protein sequences downloaded from Ensembl and output the proteinID and function in two columns
# usage "pareProteinFunction.py [input file]"

import sys
import re

# getting input path
try:
    inputFile = sys.argv[1]
    outputFile = sys.argv[2]
    print "inputFile", inputFile
    print "outputFile", outputFile
except:
    print "No inputFile and/or outputFile provided"
    print "Usage: python parseProteinFunction.py [inputFile] [outputFile]"
    sys.exit(1)

inFile = open(inputFile, 'r')
outFile = open(outputFile, 'w')

for line in inFile:
    characters = line.split(',')
    print "line", line
    protID = list(characters[0]) # "ENSP00000171887"
    description = list(characters[1])
    print "characters[0]", characters[0]
    print "characters[1]", characters[1]
    print "protID[:4]", protID[:3]
    if protID[:3] == ['E', 'N', 'S']:
        writeLine = characters[0] + "\t" + characters[1] + "\n"
        outFile.write(writeLine)
        print "writeLine", writeLine 

outFile.close()

