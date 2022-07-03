#!/usr/bin/python
# Jun Hoe, Lee (2012)
# view statistical properties of a bed file
# to call: bedStat file database [-options]
# Input: single bed file, combined cat of several bed files,
# (e.g. cat file1.bed file2.bed | bedStat stdin database)
# Input data can be unsorted, script should sort data
# Output: terminal, output data should be arranged to to make eyeballing easy

from __future__ import division  # forces division to return floating values
from operator import itemgetter  # for sorting the input bedfile
import sys
import os
import subprocess
import argparse
import tempfile
import re
import copy
import operator

import numpy as np
import matplotlib.pyplot as plt


#############################################
# class to display help message if there is error in parameter
class errorDisplayParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

# parser for bedfile input, database and optional arguments
parser = errorDisplayParser(description='bedStat: Quick statistics for a bed file.')
parser.add_argument('bedfile', type=argparse.FileType('rt'),
                    help='bed file input by user')  # sys.argv[1] bedfile input
parser.add_argument('database', action='store',
                    help='genome build, e.g danRer7')  # database
parser.add_argument('-percentile', type=float, action='store',
                    help='percentile of element pairs sorted according to spacing size. Usage: -percentile [0.01 - 1.00]')
parser.add_argument('-spacingOut', type=argparse.FileType('w'),
                    help='create bed file that lists coordinates of all pairs/triplets/etc according to -spacingOutDist. Usage: -spacingOut [filename]')
parser.add_argument('-spacingOutDist', action='store', default='<10',
                    help='specifies limit of spacing distance for element pairs listed in -spacingOut (default: \'<10\'). Usage: -spacingOutdist [conditional operator + number] (accepted number: >=1)')
parser.add_argument('-overlapOut', type=argparse.FileType('w'),
                    help='create bed file that lists coordinates of all pairs/triplets according -overlapOutDist. Usage: -overlapOut [filename]')
parser.add_argument('-overlapOutDist', action='store', default='a>=1',
                    help='''specifies limit of absolute/relative overlap distance for -overlapOut (default: \'a>=1\'). 
                          Usage1: -overlapOutDist [a + conditional operator + number] (accepted number: >=1).
                          Usage2: -overlapOutDist [r + conditional operator + percent] (accepted percent: 0<x<1.0)''')
parser.add_argument('-winSize', type=int, action='store', default='1000000',
                    help='specifies size of sliding window (default: 1 Mbp). Usage: -winSize [number]')
parser.add_argument('-offset', type=int, action='store', default='100000',
                    help='number of base positions to slide the window (default: 100 kbp). Usage: -offset [number]')
parser.add_argument('-winOut', type=argparse.FileType('w'),
                    help='create bed file that lists coordinates of windows that contain number of elements of -winOutNum. Usage: -winOut [filename]')
parser.add_argument('-winOutNum', action='store', default='>=1',
                    help='specifies limit of number of elements in a window for -winOut (default: \'>=1\'). Usage: -winOutNum [conditional operator + number]')
parser.add_argument('-keepTmpFiles', action='store_true', default=False, help='keep logged files for bedCov, sliding windows and overlapSelect. Usage: -keepTmpfiles')  # keep temp files flag
parser.add_argument('-verbose', action='store_true', default=False, help='turns on verbose mode. Usage: -verbose')  # verbose flag
parser.add_argument('-noHist', action='store_false', default=True, help='turns off the displaying of histograms. Usage: -noHist')  # noHist flag


args = parser.parse_args()

############################################
# list of functions to be called in this program

# function to turn on verbose mode
if args.verbose:  # function turned on if verbose called
    def verboseprint(arg, *args):  # *args to output any/all given arguments
        print "-v ", arg, ":", args
else:
    verboseprint = lambda *a: None  # do-nothing function


# flag to keep tempFiles for bedCov, sliding windows and overlapSelect
keepTemp = True  # default delete option for tempFiles
if args.keepTmpFiles:
    keepTemp = False


# function to calculate percentile: # called to calculate q1, med, q3 or manual percentile
def percentile(sample_list, percent):
    """
    @parameter sample_list: list of values, must be sorted
    @parameter percent: specified percentage, e.g. q1, med, q3
    """
    verboseprint('sorted sample list', sample_list)
    percent_pos = int(round(percent * len(sample_list) + 0.5))  # position of q1
    percent_value = sample_list[percent_pos - 1]  # value of Q1
    verboseprint('percentile, position of percentile, number of samples', percent_value, percent_pos, len(sample_list))
    return percent_value


# function to calculate average
def average(sample_list):
    """
    @parameter sample_list: list of values
    """
    avg = sum(sample_list) / len(sample_list)
    verboseprint('average (sum/number), sum of element sizes, number of elements',
                  avg, sum(sample_list), len(sample_list))
    return avg


# function to create histogram subplot (for spacing and overlap)
def histogram(fig, y, N, barColour, ticks_x, label_x):
    """
    @parameter fig: subplot number
    @parameter y: dataset (in list)
    @parameter N: number of samples, len(y)
    @parameter barColour: colour of bars, e.g. 'r'
    @parameter ticks_x: list of items for x-axis
    @parameter label_x: label for x-axis
    """
    plt.subplot(2, 2, fig)  # 2 rows, 2 column, subplot/figure number (1,2,3..)
    x = np.arange(1, N + 1)
    width = 1.0
    plt.bar(x, y, width, color=barColour)
    plt.xticks(x + (width / 2), ticks_x)
    plt.xlabel(label_x)
    plt.ylabel('Number of elements')


# function to pass an operator as a parameter (for spacingOut, overlapOut)
def compare(itemInList, function, outNumber):
    """
    @parameter itemInList: item in space_list or winList that will be iterated
    @parameter function: operator, given by user ('<', '<=', '>', '>=')
    @parameter outNumber: value given by user for spacingOut, overlapOut, winOut
    """
    mappings = {'<': operator.lt, '<=': operator.le,
               '>': operator.gt, '>=': operator.ge}  # map each operator to a function
    verboseprint('compare function: Item in list, operator, limit value', itemInList, function, outNumber)
    return mappings[function](itemInList,outNumber)


# function to output consecutive items within a range (cluster), (for spacingOut, overlapOut)
def cluster(given_list, inLimit, number):
    """
    @given_list: list to be iterated, i.e. space_list or winList
    @inLimit: operator for items in list within value, pass to function 'compare'
    @number: value given by user for spacingOut, overlapOut, winOut to be passed to function 'compare'
    """
    inCluster = False  # flag to check whether an item is in cluster
    start = end = None # start = position (index) for chromStart; end = position (index) for chromEnd
    for i in range(1, len(given_list)):
        verboseprint('index, space, start, end', i, given_list[i-1], start, end)
        if not compare(given_list[i], inLimit, number) or given_list[i] == '-':  # if compare function returns False or item is == '-':
            if inCluster:
                end = i - 1
                verboseprint('Break cluster: index, space, start, end', i, given_list[i-1], start, end)
                yield start, end
                inCluster = False  # reset inCluster flag to False
        else:  # if compare function returns True
            if not inCluster:  # if inCluster flag is False
                start = i - 1
                inCluster = True  # set inCluster flag to True 
                verboseprint('Start cluster: index, space, start', i, given_list[i-1], start) 


# function to equalize '<=/>=' into '</>' for -spacingOutDist and -overlapOutDist
def condition(conditional, number):  
    """
    @parameter conditional: Conditional operator in -spacingOutDist and -overlapOutDist
    @parameter number: Number/value in -spacingOutDist and -overlapOutDist
    """
    if conditional == '<=':
        number = number + 1  # converts >= to >, by subtracting 1
    elif conditional == '>=':
        number = number - 1  # converts >= to >, by adding 1
    return number


# function to compare between two element/window pairs (for use in clustering/windows)
# input parameters are list of (number of element, chrom, chromStart, chromEnd)
def isElementPairOverlapping(winExam, winCompare): 
    """
    @parameter winExam: the window being examined
    @parameter winCompare: the window being used as comparison
    """
    if winExam[1] == winCompare[1]: # second item in the list of tuples (number of element, chrom, chromStart, chromEnd)
        if winExam[3] <= winCompare[2] or winExam[2] >= winCompare[3]:  # NOTE that in bed 'end' is the base downstream of the last base in the element --> therefore >= and <= comparisons
            return False  # False == no overlap found
        else:
            return True  # True == overlap found
    else:
        return False  # False == no overlap found


# function to compare whether an element overlaps with another element in a list (for use in clustering/windows)
def isElementOverlappingWithList(sample_list, selected_list):
    """
    @parameter sample_list  : the list of sample windows to be checked
    @parameter selected_list: list of windows that have been chosen and being used as comparison
    """
    for i in sample_list:
        overlapCheckCounter = 0  # check for overlap occurences, resets for each item in sample_list
        if len(selected_list) >= 5: # break after 5 items have been selected
            break
        for j in selected_list:
            checkOverlap = isElementPairOverlapping(i,j)  # call function isElementPairOverlapping
            if checkOverlap:  # if checkOverlap returns True (i.e. overlap found)
                overlapCheckCounter += 1
        if overlapCheckCounter == 0:  # if checkOverlap returns False against every item in selected_list
            selected_list.append(i)

#############################################
### checks for user input bedfile and values

# check if valid database name is provided by user
chromSizePath = '/genome/gbdb-HL/'
chromSizeFile = chromSizePath + args.database + '/chrom.sizes'  # join items to form path to /chrom.sizes
if not os.path.exists(chromSizeFile):  # if no valid path to database is found
    print 'Incorrect or invalid database (genome build) provided.'
    print 'usage: python bedStat.py bedFile database [options]'
    sys.exit(1)

# check arguments given for -winSize, -offSet
if args.winSize < 1:
   print 'Incorrect value provided for -winSize.'
   sys.exit(1)
if args.offset < 1:
   print 'Incorrect value provided for -offSet.'
   sys.exit(1)

############################################
### process input bedfile through sorting, calculating element sizes and spacing between elements

# reading input bed file from user and convert each element (line) into a tuple within a master list, e.g. [(chr1, 100, 200), (chr1, 300, 400)..]
# not using the kentTool bedSort as it only sorts by chrom, chromStart and can't distinguish for elements with identical chromStart but different chromEnd
inbedSort = [] # list to store the element tuples
try:
    for lineNum, line in enumerate(args.bedfile, start=1):  # lineNum used if error is raised
        if not line.strip():  # skip blank lines
            continue
        elif line.strip():
            bedline = line.split('\t')
            element = (bedline[0], int(bedline[1]), int(bedline[2]))  # have to combine them into a tuple first as list.append takes only 1 argument
            inbedSort.append(element) # puts the element tuple into the outbedSort.list
except:
    print "Error: Cannot read line", lineNum, ". Please check bedfile."  # raise exception if there are problems with data in the bedfile
    print line
    sys.exit(1)

sortedbedSort = sorted(inbedSort, key=itemgetter(0,1,2)) # sort tuples first by chrom, then chromStart, finally chromEnd

# create temporary bedfile outbedSort to store the output of chrom, chromStart, chromEnd
# outbedSort will be used subsequently as the input file for bedCov and overlapSelect
outbedSort = tempfile.NamedTemporaryFile(prefix='bedSort', delete=keepTemp)  

# now unpack the tuples and sort them into lists of chrom, chromStart, chromEnd
chrom_list = []  # list of chrom number/groups
chromStart_list = []  # list of chromStart positions
chromEnd_list = []  # list of chromEnd positions

for chrom, chromStart, chromEnd in sortedbedSort:  # directly reference the three items in each tuple of sortedbedSort
    if not line.strip():  # skip blank lines
        continue
    else:
        chrom_list.append(chrom)
        chromStart_list.append(chromStart)
        chromEnd_list.append(chromEnd)
        outbedSort.write('%s\t%s\t%s\n' % (chrom, chromStart, chromEnd)) #  write to outbedSort

outbedSort.seek(0)  # resets cursor position for reading by bedCov

# calculate eleSize spacing & overlap (negative spacing)
eleSize_list = []  # list of element sizes
space_list = ['-']  # list of spaces and overlaps (negative spaces), starts with '-' to ensure len(chrom_list[i]) = chrom(Start_list[i])

eleSize_list.append(chromEnd_list[0] - chromStart_list[0]) # adding the first item to eleSize_list as space_list calculation starts from the second element
for i in range(1, len(chromStart_list)):
    if chrom_list[i] == chrom_list[i - 1]:  # within the same chromosome group
        if chromEnd_list[i] <= chromEnd_list[i - 1]:  # if subsequent element is smaller than previous element
            space = chromStart_list[i] - chromEnd_list[i]  # whereby overlap = ele_Size
        else:
            space = chromStart_list[i] - chromEnd_list[i - 1]
        space_list.append(space)
    else:
        space_list.append('-')  # if different chromo group, append '-' to ensure chrom_list[i] = chromStart_list[i] = space_list[i]
    eleSize_list.append(chromEnd_list[i] - chromStart_list[i])
space_list.append('-')  # add extra dummy item to space_list
space_len = len(space_list)

if args.verbose:
    print '-v i(index), chrom, chromStart, chromEnd, eleSize, spacing(chromStart(i+1) - chromEnd(i):'
    for i in range(len(chrom_list)):
        print i, chrom_list[i], chromStart_list[i], chromEnd_list[i], eleSize_list[i], space_list[i]

# create list of element pairs that have spaces(spaceAbles), i.e. excluding '-' from different chromosomes
spaceAbles = space_len - space_list.count('-') 
verboseprint('check: len(chrom_list) == len(chromStart_list) == len(chromEnd_list)', len(chrom_list), len(chromStart_list), len(chromEnd_list))
if args.verbose:
    print "-v * finished reading input bedfile. Created lists of chrom, chromStart, chromEnd, eleSize, spacing\n"

#############################################################################
###A. size information
if args.verbose:
    print "-v A. Statistics/size information"

sortEleSize_list = copy.deepcopy(eleSize_list)  # make a deep copy of eleSize_list
sortEleSize_list.sort()  # sorts lists, original eleSize_list position of items is unchanged
bedStatq1 = percentile(sortEleSize_list, 0.25)  # calling percentile function to calculate q1
bedStatmed = percentile(sortEleSize_list, 0.50)
bedStatq3 = percentile(sortEleSize_list, 0.75)
bedMean = average(eleSize_list)
print 'Statistics for size of elements (bp):\t\
       q1: %-12smed: %-12sq3: %-12smean: %-12.3f \
       min: %-12s max: %-12s numElements: %d' % (bedStatq1, bedStatmed, bedStatq3, bedMean, 
                                                 min(eleSize_list), max(eleSize_list), len(eleSize_list))

# if user invokes -percentile option
if args.percentile:
    if args.percentile > 0 and args.percentile < 1:
        user_percentile = percentile(sortEeleSize_list, args.percentile)  # call percentile function
        print "User percentile(%s):\t%s" % (args.percentile, user_percentile)
    else:
        print "Incorrect value provided for percentile. Accepted values: <0 < x <1.0"

if args.verbose:
    print '-v * finished section on statistics of input file\n'

##############################################################################
###B. genomic overlap.  How many bases in the genome are covered by at least one element
if args.verbose:
    print "-v B. bedCov/genomic overlap"

# create temporary file to log output of bedCov
outbedCov = tempfile.NamedTemporaryFile(prefix='bedCov', delete=keepTemp)  # default of keepTemp is True unless args.keepTmpFiles is called

try:
    bedCovIntro = "\ngenomic overlap (bedCov): "
    procBedCov = subprocess.Popen("bash -c \"bedCov " + args.database + " " + outbedSort.name
                                  + " 2>&1 | tee " + outbedCov.name + " \"", stdout = subprocess.PIPE, shell=True)
    bedCovOutput = procBedCov.communicate()
    showOutput = bedCovIntro + bedCovOutput[0]
    print showOutput
except:
    print 'bedCov cannot read bed file'

if args.verbose:
    print '-v * finished section on bedCov\n'

########################################################################################
###C. spacing between elements
print "SPACING - distance (base pairs) between two subsequent elements (within the same chromosome/scaffold)"
if args.verbose:
    print "-v C. Spacing section starts"

# i. filtering items in space_list < y, y = 0, 10, 100, 1000
# create lists to store items
spaceL0 = []  # spacing: <0 , i.e. overlaps
space0 = []  # spacing = 0
spaceL10 = []  # spacing: 0 < x <=10
spaceL100 = []  # spacing: 10 < x <= 100
spaceL1k = []  # spacing: 100 < x <= 1000
spaceL10k = []  # spacing: 1000 < x <= 10000
spaceM10k = []  # spacing: > 10,000
spaceM100k = []  # spacing: > 100,000
spaceM1mb = []  # spacing:> 1,000,000
spaceM10mb = []  # spacing: > 10,000,000

for i, elem in zip(range(space_len), space_list):
    if elem < 0:  # filtering elements with spacing == 0
        spaceL0.append(elem)
    if elem == 0:
        space0.append(elem)
    if elem > 0 and elem <= 10:  # filtering element sizes <= 10
        spaceL10.append(elem)
    if elem > 10 and elem <= 100:
        spaceL100.append(elem)
    if elem > 100 and elem <= 1000:
        spaceL1k.append(elem)
    if elem > 1000 and elem <= 10000:
        spaceL10k.append(elem)
    if elem > 10000 and elem <= 100000:  # filtering element sizes > 10k, but != '-' to avoid blanks
        spaceM10k.append(elem)
    if elem > 100000 and elem <= 1000000:
        spaceM100k.append(elem)
    if elem > 1000000 and elem <= 10000000:
        spaceM1mb.append(elem)
    if elem > 10000000 and elem != '-':
        spaceM10mb.append(elem)

verboseprint('spaceL0 (overlaps), len(spaceL0)', spaceL0, len(spaceL0))
verboseprint('space0 (spacing==0), len(space0)', space0, len(space0))
verboseprint('spaceL10 (0 < x <=10), len(spaceL10)', spaceL10, len(spaceL10))
verboseprint('spacel100 (10 < x <= 100) , len(spaceL100)', spaceL100, len(spaceL100))
verboseprint('spaceL1 k (100 < x <= 1000), len(spaceL1k)', spaceL1k, len(spaceL1k))
verboseprint('spaceL10 k (1000 < x <= 10000), len(spaceL10k)', spaceL10k, len(spaceL10k))
verboseprint('spaceM10 k (> 10,000), len(spaceM10k)', spaceM10k, len(spaceM10k))
verboseprint('spaceM100 k (> 100,000), len(spaceM100k)', spaceM100k, len(spaceM100k))
verboseprint('spaceM1 M (> 1,000,000), len(spaceM1mb)', spaceM1mb, len(spaceM1mb))
verboseprint('spaceM10mb (> 10,000,000), len(spaceM1mb)', spaceM10mb, len(spaceM1mb))
verboseprint('Check: total items with spaces == spaceAbles', 
             len(spaceL0) + len(space0) + len(spaceL10) + len(spaceL100) + len(spaceL1k) + len(spaceL10k) + len(spaceM10k) + len(spaceM100k) + len(spaceM1mb) + len(spaceM10mb),
             spaceAbles)

print "E.g.\"chr1	100	200, chr1	350	400\". Spacing: 350 - 200 = 150 bp"
print "Absolute spacing - number of elements pairs that possess the displayed spacing size (base pairs)."
print "Percentage - percentage of element pairs with the displayed spacing size out of the total number of element pairs." 
print "Spacing (absolute spacing, percentage)"
print "\t<0 bp (overlaps): %15s, %10.3f%%" % (len(spaceL0), len(spaceL0) * 100 / spaceAbles)  # absolute space, relative space
print "\t0 bp            : %15s, %10.3f%%" % (len(space0), len(space0) * 100 / spaceAbles)
print "\t<=10 bp         : %15s, %10.3f%%" % (len(spaceL10), len(spaceL10) * 100 / spaceAbles)
print "\t<=100 bp        : %15s, %10.3f%%" % (len(spaceL100), len(spaceL100) * 100 / spaceAbles)
print "\t<=1 kbp         : %15s, %10.3f%%" % (len(spaceL1k), len(spaceL1k) * 100 / spaceAbles)
print "\t<=10 kbp        : %15s, %10.3f%%" % (len(spaceL10k), len(spaceL10k) * 100 / spaceAbles) 
print "\t>10 kbp         : %15s, %10.3f%%" % (len(spaceM10k), len(spaceM10k) * 100 / spaceAbles)
print "\t>100 kbp        : %15s, %10.3f%%" % (len(spaceM100k), len(spaceM100k) * 100 / spaceAbles)
print "\t>1 Mbp          : %15s, %10.3f%%" % (len(spaceM1mb), len(spaceM1mb) * 100 / spaceAbles)
print "\t>10 Mbp         : %15s, %10.3f%%" % (len(spaceM10mb), len(spaceM10mb) * 100 / spaceAbles)

# ii. create spacing histogram, spaceHist by calling histogram function
spaceHist_y = [len(spaceL0), len(space0), len(spaceL10), len(spaceL100), len(spaceL1k), 
               len(spaceL10k), len(spaceM10k), len(spaceM100k), len(spaceM1mb), len(spaceM10mb)]
spaceHist_x = ['<0', '0', '<=10', '<=100', '<=1 k', '<=10 k', '>10 k', '>100 k', '>1 M', '>10 M']
# call histogram function: (figure 1, spacing data for y-axis, 
#                           number of samples, red colour for bars, bar labels for x-axis, label for x-axis)
spaceHist = histogram(1, spaceHist_y, len(spaceHist_y), 'r', spaceHist_x, 'Size(spacing, bp)')

# iii. create outfile if user invokes '-spacingOut [filename]'
regexspacingOutDist = re.compile(r'([><=]+)([0-9]+)')
spaceMatch = regexspacingOutDist.search(args.spacingOutDist)
if args.spacingOut:  # if '-spacingOut' option for output file is called
    if spaceMatch:  # if match exists
        spaceCondition = spaceMatch.group(1)  # conditional operator, e.g. '>'
        spaceNumber = spaceMatch.group(2)  # number, e.g. '50'
        spaceNumber = int(spaceNumber)  # converts spaceNumber to integer 
        if spaceNumber >= 1:  # checks minimum value for spaceNumber
            # distLoop = counter for internal loop, for cases where the subsequent element is also < spaceNumber
            # e.g. space [i] < spaceNumber; space[i+1] < spaceNumber, ...
            # if subsequent elements < spaceNumber, distLoop prevents "for i in range(space_len)" from iterating over the line again.
            if spaceCondition == '<' or spaceCondition == '<=':
                if args.verbose:
                    print '''-v spacingOut. How the algorithm works: (e.g. <100)
                             E.g. chr1	100	200	-
                                  chr1	250	300	50(spacing)
                                  chr1	375	400     75
                                  chr2	800	825	-
                             Entries considered: line 2, line3, line4 (loop breaks)
                             File output: chr1	100	400'''
                for start, end in cluster(space_list, spaceCondition, spaceNumber): # call cluster function
                    args.spacingOut.write('%s\t%s\t%s\n' % (chrom_list[start],
                                                            chromStart_list[start], chromEnd_list[end]))
		    verboseprint('spacingOut output: index, chrom, chromStart, chromEnd, space (of last element in cluster)',
                                 i, chrom_list[start], chromStart_list[start], chromEnd_list[end], space_list[end])
                print "Note: Subsequent elements that fulfil spacing criteria are combined into a single bed entry, e.g: \'<50\'"
                print "\"chr1 100-200 bp, chr 1 210-250, chr1 280-300\" -> \"chr1	100	300\""
                args.spacingOut.close()
            elif spaceCondition == '>' or spaceCondition == '>=':
                if args.verbose:
                    print '''-v spacingOut. How the algorithm works: (e.g. > 100)
                                E.g. chr1	100	200	-
                                     chr1	350	400	150
                                     chr2       450	550	-
                             Entries considered: line 2, line 3 (loop breaks)
                             File output: chr1	100	400'''
                for i, item in zip(range(space_len), space_list):
                    if compare(item, spaceCondition, spaceNumber) and item != '-':  # call compare function, if returns True
                        args.spacingOut.write('%s\t%s\t%s\n' % (chrom_list[i], chromStart_list[i-1], chromEnd_list[i]))
                        verboseprint('spacingOut output: index, chrom, chromStart, chromEnd, space',
                                     chrom_list[i], chromStart_list[i - 1], chromEnd_list[i], space_list[i])
                args.spacingOut.close()
    else:  # args.spacingOutDist != spaceMatch.group(0)
        print "Incorrect input values given for -spacingOutDist. Output file process aborted."

if args.verbose:
    print '-v *finished section on spacing'

###############################################################################
###D. overlap
print "\nOVERLAP - two or more elements that possess overlapping positions (e.g. chr1:100-200, chr1:150-300)"
if args.verbose:
    print "D. Overlap section"

# i. calculate absolute overlap, where overlap = -space(negative space)
# create lists to store values for absolute overlaps
notOverlap0 = []  # not overlaps, i.e. has spacing >=0
absOverlapL10 = []  # spacing: -10 <= x < 0,
absOverlap10_100 = []  # spacing: -100 <= x < -10
absOverlap100_1k = []  # spacing: -1000 <= x < -100
absOverlapM1k = []  # spacing: > -1000

# create lists to store values for relative overlaps (e.g. element 1:100-200, element 2:150-300, relative overlap=50/100)
relOverlap_list = []
relOverlap10 = [] # relOverlap: <= 10
relOverlap10_30 = []  # relOverlap: 10 < x <= 30
relOverlap30_50 = []  # relOverlap: 30 < x <= 50
relOverlap50_80 = []  # relOverlap: 50 < x <= 80
relOverlap80_95 = []  # relOverlap: 80 < x <= 95
relOverlap95_100 = []  # relOverlap: 95 < x <= 100
if args.verbose:
    relOverlap100 = []  # relOverlap: > 100

# get values for absolute overlaps and create the relative overlap list
# elements with 0 spacing is not counted as overlaps, but as spacing = 0
if args.verbose:
    print '''Calculating relOverlap:
            E.g. chr1	100	200
                 chr1	150	400
                spacing = 150-200 , smallerElement = (200-100)
                relOvelap = positive(-50) x 100/smallerElement'''
for i, elem in zip(range(space_len), space_list):
    if elem >= 0 and elem != '-':  # filtering spaceAble elements that don't overlap,
        notOverlap0.append(elem)
    if elem < 0 and elem >= -10:  # filtering elements with overlap 0>x>=10
        absOverlapL10.append(elem)
    if elem < -10 and elem >= -100:  # filtering element overlaps < 100
        absOverlap10_100.append(elem)
    if elem < -100 and elem >= -1000: # filtering element overlaps < 100
        absOverlap100_1k.append(elem)
    if elem < -1000:
        absOverlapM1k.append(elem)
    if elem < 0:  # making list of elements that overlap to calculate relative overlap
        smallerElement = int(min(eleSize_list[i], eleSize_list[i - 1]))  # smaller size of two elements
        relOverlap = (abs(elem) * float(100)) / smallerElement
        relOverlap_list.append(relOverlap)
        verboseprint('i, relOvelap, space, eleSize_list[i-1], eleSize_list[i], smallerElement',
                      i, relOverlap, elem, eleSize_list[i-1], eleSize_list[i], smallerElement)
    if elem >= 0: # filtering elements that have spaces or just '-'
        relOverlap_list.append('-')  # ensures that len(relOverlaplist) == len(chrom_list) for chrom positions

verboseprint('notOverlap0 (spacing), len(notOverlap0)', notOverlap0, len(notOverlap0))
verboseprint('absOverlapL10 (-10 <= x < 0), len(absOverlapL10)', absOverlapL10, len(absOverlapL10))
verboseprint('absOverlap10_100 (-100 <= x < -10), len(absOverlap10_100)', absOverlap10_100, len(absOverlap10_100))
verboseprint('absOverlap100_1k (-1000 <= x < -100), len(absOverlap100_1k)', absOverlap100_1k, len(absOverlap100_1k))
verboseprint('absOverlapM1k (> -1000), len(absOverlapM1k)', absOverlapM1k, len(absOverlapM1k))
verboseprint('Check: total items of absolute overlaps = spaceAbles - notOverlap0',
             len(absOverlapL10) + len(absOverlap10_100) + len(absOverlap100_1k) + len(absOverlapM1k), spaceAbles, len(notOverlap0))

for i, elem in zip(range(len(relOverlap_list)), relOverlap_list):  # get values for relative overlap
    if elem <= 10:
        relOverlap10.append(elem)
    elif elem > 10 and elem <= 30:
        relOverlap10_30.append(elem)
    elif elem > 30 and elem <= 50:
        relOverlap30_50.append(elem)
    elif elem > 50 and elem <= 80:
        relOverlap50_80.append(elem)
    elif elem > 80 and elem <= 95:
        relOverlap80_95.append(elem)
    elif elem > 95 and elem <= 100:
        relOverlap95_100.append(elem)
    if args.verbose:
        if elem > 100 and elem != '-':
            relOverlap100.append(elem)

verboseprint('relOverlap10, len(relOverlap10)', relOverlap10, len(relOverlap10))
verboseprint('relOverlap10_30, len(relOverlap10_30)', relOverlap10_30, len(relOverlap10_30))
verboseprint('relOverlap30_50, len(relOverlap30_50)', relOverlap30_50, len(relOverlap30_50))
verboseprint('relOverlap50_80, len(relOverlap50_80)', relOverlap50_80, len(relOverlap50_80))
verboseprint('relOverlap80_95, len(relOverlap80_95)', relOverlap80_95, len(relOverlap80_95))
verboseprint('relOverlap95_100, len(relOverlap95_100)', relOverlap95_100, len(relOverlap95_100))
if args.verbose:
    print "-v Check any non-intended relative overlap > 100%% %17s, %10.3f%%" % (len(relOverlap100), len(relOverlap100) * 100 / spaceAbles)
verboseprint('Check: total items of relative overlaps (excluding spaceL0) = spaceAbles - notOverlap0', 
             len(relOverlap10) + len(relOverlap10_30) + len(relOverlap30_50) + len(relOverlap50_80) + len(relOverlap80_95) + len(relOverlap95_100), spaceAbles, len(notOverlap0))
verboseprint('check: total items absOverlap == total relOverlap',
             len(absOverlapL10) + len(absOverlap10_100) + len(absOverlap100_1k) + len(absOverlapM1k),
             len(relOverlap10) + len(relOverlap10_30) + len(relOverlap30_50) + len(relOverlap50_80) + len(relOverlap80_95) + len(relOverlap95_100))
verboseprint('check: spaceAbles = total items with spaces + total item with overlaps', spaceAbles,
             len(space0) + len(spaceL10) + len(spaceL100) + len(spaceL1k) + len(spaceL10k) + len(spaceM10k) + len(spaceM100k) + len(spaceM1mb) + len(spaceM10mb),
             len(absOverlapL10) + len(absOverlap10_100) + len(absOverlap100_1k) + len(absOverlapM1k))
verboseprint('check: number of \'-\' in relOverlap_list (non-overlaps) = space_len - total items with relOverlaps', relOverlap_list.count('-'), space_len,
             len(relOverlap10) + len(relOverlap10_30) + len(relOverlap30_50) + len(relOverlap50_80) + len(relOverlap80_95) + len(relOverlap95_100))

if any(elem != '-' for elem in relOverlap_list):  # if overlapping elements are found, in which items other than '-' are found in relOverlap_list
    print "1. Absolute overlap - number of two (or more) elements that possess the following overlapping positions."
    print "Percentage - percentage of elements with the following overlap size out of the total number of element pairs."
    print "Absolute overlap (absolute overlap, percentage)"  # printing output for absolute overlap
    print "\t<=0 bp (spaces): %15s, %10.3f%%" % (len(notOverlap0), len(notOverlap0) * 100 / spaceAbles)  # absolute overlap, percentage
    print "\t<=10 bp        : %15s, %10.3f%%" % (len(absOverlapL10), len(absOverlapL10) * 100 / spaceAbles) 
    print "\t10< x<=100 bp  : %15s, %10.3f%%" % (len(absOverlap10_100), len(absOverlap10_100) * 100 / spaceAbles)
    print "\t100< x<=1000 bp: %15s, %10.3f%%" % (len(absOverlap100_1k), len(absOverlap100_1k) * 100 / spaceAbles)
    print "\t>1000 bp       : %15s, %10.3f%%" % (len(absOverlapM1k), len(absOverlapM1k) * 100 / spaceAbles)

    print "\n2. Relative overlap (Overlap/smaller element size)"  # printing output for relative overlap
    print "\tE.g. chr1 	100	200"
    print "\t     chr1	150	400"
    print "\tRelative overlap: (200-150)/100 = 50%"
    print "\t<=0 bp (non-overlapping)%%: %15s, %10.3f%%" % (len(notOverlap0), len(notOverlap0) * 100 / spaceAbles)  # double %% as escape character for %
    print "\t0  < x <= 10%%   : %15s, %10.3f%%" % (len(relOverlap10), len(relOverlap10) * 100 / spaceAbles)  # double %% as escape character for %
    print "\t10 < x <= 30%%   : %15s, %10.3f%% " % (len(relOverlap10_30), len(relOverlap10_30) * 100 / spaceAbles)
    print "\t30 < x <= 50%%   : %15s, %10.3f%%" % (len(relOverlap30_50), len(relOverlap30_50) * 100 / spaceAbles)
    print "\t50 < x <= 80%%   : %15s, %10.3f%%" % (len(relOverlap50_80), len(relOverlap50_80) * 100 /spaceAbles)
    print "\t80 < x <= 95%%   : %15s, %10.3f%%" % (len(relOverlap80_95), len(relOverlap80_95) * 100 / spaceAbles)
    print "\t95 < x <=100%%   : %15s, %10.3f%%" % (len(relOverlap95_100), len(relOverlap95_100) * 100 / spaceAbles)
else:  # no overlapping elements found
    print "No overlap found between any element pairs."

# iii. create histograms
# histogram of absolute overlap values
absoverlapHist_x = ['<=0', '<=10', '<=100', '<=1000', '>1000']  # x-axis, category of element sizes
absoverlapHist_y = [len(notOverlap0), len(absOverlapL10), len(absOverlap10_100), len(absOverlap100_1k), len(absOverlapM1k)]  # y-axis, number of elements
# call histogram function: (figure 2, abs_overlap data for y-axis, number of samples,
#                           blue colour for bars, bar labels for x-axis, label for x-axis)
absoverlap_hist = histogram(2, absoverlapHist_y, len(absoverlapHist_y), 'b', absoverlapHist_x, 'Size(absolute overlap, bp)')

# histogram of relative overlap values
reloverlapHist_x = ['<=10%', '<=30%', '<=50%', '<=80%', '<=95%', '<=100%']
reloverlapHist_y = [len(relOverlap10), len(relOverlap10_30), len(relOverlap30_50),len(relOverlap50_80),
                    len(relOverlap80_95), len(relOverlap95_100)]
# call histogram: (figure 3, rel_overlap data for y-axis, number of samples,
#                  green colour for bars, bar labels for x-axis, label for x-axis)
reloverlap_hist = histogram(3, reloverlapHist_y, len(reloverlapHist_y), 'g', reloverlapHist_x, 'Size(relative overlap, bp)')

# iv. create outfile, if user invokes '-overlapOut'
regexoverlapOutDist = re.compile(r'([ar])([><=]+)([\d\W\d]+)')
overlapMatch = regexoverlapOutDist.search(args.overlapOutDist)
if args.overlapOut:  # if user calls optional output 'overlapOut'
    if overlapMatch:  # if match exists, i.e. overlapMatch returns 'True'
        overlapType = overlapMatch.group(1)  # get type of overlap, absolute or relative
        overlapCondition = overlapMatch.group(2)  # get conditinal operator
        overlapNumber = overlapMatch.group(3)  # get value of overlap
        if overlapType == 'a':  # chunk code for absolute relative
           try:
               overlapNumber = int(overlapNumber)  # try converting overlapNumber to integer 
           except:
               print "Incorrect input value for parameter -overlapOutDist. Output file process aborted."
           if isinstance(overlapNumber, int):  # if overlapNumber is successfully converted to integer
               overlapNumber = condition(overlapCondition, overlapNumber) # call conditional function
               if overlapCondition == '<' or overlapCondition == '<=':
                   if args.verbose:
                       print '''-v overlapOut. How the algorithm works: (e.g. < 100)
                                E.g. chr1 100	300	-
                                     chr1 250	400	-50
				     chr2 450	550	-
                                Entries considered: line 2, line 3 (loop breaks)
                                File output: chr1	100	400'''
                   for i, item in zip(range(space_len), space_list):
                      if compare(item, '>', -overlapNumber) and compare(item, '<', 0):  # call compare function, if return True
                          args.overlapOut.write('%s\t%s\t%s\n' % (chrom_list[i], chromStart_list[i-1], chromEnd_list[i]))
                          verboseprint('overlapOut output: index, chrom, chromStart, chromEnd, space',
                                       chrom_list[i], chromStart_list[i-1], chromEnd_list[i], space_list[i]) 
                   args.overlapOut.close()
               elif overlapCondition == '>' or overlapCondition == '>=':
                   if args.verbose:
                      print '''-v overlapOut. How the algorithm works: (e.g. a > 20)
                                  E.g. chr1	100	200	-
		                       chr1	150	200	-50(spacing)
                                       chr1	175	300     -25
                                       chr2	800	825	-
                                  Entries considered: line 2, line3, line4 (loop breaks)
                                  File output: chr1	100	400'''
                   for start, end in cluster(space_list, '<', -overlapNumber): # call cluster function
                       args.overlapOut.write('%s\t%s\t%s\n' % (chrom_list[start],
                                             chromStart_list[start], chromEnd_list[end]))
		       verboseprint('overlapOut output: index, chrom, chromStart, chromEnd, space (of last element in cluster)',
                                     i, chrom_list[start], chromStart_list[start], chromEnd_list[end], space_list[end])
                   args.overlapOut.close()
                   print "Note: Subsequent elements that fulfil overlap criteria are combined into a single bed entry, e.g:  <50"
                   print "\"chr1 100-200 bp, chr 1 100-250, chr1 200-300\" -> \"chr1	100	300\""
        elif overlapType == 'r' :  # chunck code for relative overlap
            try:
                overlapNumber = float(overlapNumber)  # try converting overlapNumber to float
            except:
                print "Incorrect input value for parameter -overlapOutDist. Output file process aborted." 
            if isinstance(overlapNumber, float):  # check whether overlapNumber is successfully converted to float
                if overlapNumber > 0 and overlapNumber <= 1:  # checks whether overlapNumber is within range of 0-1.0
                    overlapNumber = overlapNumber * 100  # relOverlap_list is in percentage (0-100%); overlapNumber is in decimal(0.0-1.0)
                    if args.verbose:
                        print '''-v overlapOut. How the algorithm works: (e.g. r < 0.50)
                                  				Spacing	Size	relOverlap
                                  E.g. chr1	100	300	-	200 	-
                                       chr1	250	400	-50	150	33.33
                                       chr2 	450	550	-	100	-
                                  Entries considered: line 2, line 3 (loop breaks)
                                  File output: chr1	100	400'''
                    if overlapCondition == '<' or overlapCondition == '<=':
                        print "mali", overlapCondition
                        for i, item in zip(range(len(relOverlap_list)), relOverlap_list):
                            if compare(item, overlapCondition, overlapNumber) and item != '-':  # call compare function, if return True
                                verboseprint('index, element considered:, chrom, chromStart, chromEnd, overlap, relOverlap',
                                             i, chrom_list[i], chromStart_list[i], chromEnd_list[i], space_list[i], relOverlap_list[i])
                                args.overlapOut.write('%s\t%s\t%s\n' % (chrom_list[i], chromStart_list[i-1], chromEnd_list[i]))
                                verboseprint('overlapOut output: index, chrom, chromStart, chromEnd, relOverlap',
                                             chrom_list[i], chromStart_list[i-1], chromEnd_list[i], space_list[i])
                    elif overlapCondition == '>' or overlapCondition == '>=': 
                        for start, end in cluster(relOverlap_list, overlapCondition, overlapNumber): # call cluster function
                            args.overlapOut.write('%s\t%s\t%s\n' % (chrom_list[start],
                                                  chromStart_list[start], chromEnd_list[end]))
		            verboseprint('overlapOut output: index, chrom, chromStart, chromEnd, space (of last element in cluster)',
                                         i, chrom_list[start], chromStart_list[start], chromEnd_list[end], space_list[end])
                    args.overlapOut.close()
    else:
        print "Incorrect input value for parameter -overlapOutDist. Output file process aborted."

if args.verbose:
    print '-v *finished section on overlap'

########################################################################################
###E. clustering/sliding window

print "\nCLUSTERING WINDOW - examining number of elements in each window"

# read /genome/gbdb-HL/[args.database]/chrom.sizes to get length of chromosomes, args.database = genome of species
# open file chrom.sizes
chromData = open(chromSizeFile, 'r')  # chromData = chrom.size file in database to be read
chromNameSize = {}  # create dictionary to store chromName and chromSize e.g. {'chr1' : 12345...}
for line in chromData:
    if not line.strip():  # skip blank lines
        continue
    else:
        chromLine = line.split('\t')  # reads line and split by tab, e.g. "chr1 81811"
        chromName = chromLine[0]
        chromSize = chromLine[1]
        chromNameSize[chromName] = int(chromSize)  # converts chromName and chromSize into a dictionary {'chromName': chromSize}
        
# Use a while loop to send output (coordinates) to temporary file, winCompilebed
winCompilebed = tempfile.NamedTemporaryFile(prefix='win', delete=keepTemp)

if args.verbose:
    print "-v Creating sliding windows"

for key,value in chromNameSize.items():
    # set chrom_len as the limit, based on the chromSize (value) of each key (chromName) in the chromNameSize dictionary
    chrom_len = chromNameSize[key]
    winStart = 0  # start position of window, resets to 0 for each chrom
    winEnd = winStart + args.winSize  # end position of window
    if args.winSize > chrom_len:  # if window size > chromosome size, skip chromosome
        verboseprint('skipped chromosome/scaffolds: winSize, chrom_size)', key, args.winSize, value)    
        continue
    while winEnd <= chrom_len:  # while window size < chromosome size
        winCompilebed.write('%s\t%s\t%s\n' % (key, winStart, winEnd))
        winStart = winStart + args.offset  # shift start window to new position
        winEnd = winEnd + args.offset  # shift end window to new position

winCompilebed.seek(0)  # reset cursor position of winCompilebed
# create temp file to store output of number of elements per window, based on the output from overlapSelect.
noOfElementCount = tempfile.NamedTemporaryFile(prefix='oSelect', delete=keepTemp)

# call overlapSelect. Include options -selectFmt and -inFmt=bed because the select and input files do not have .bed suffix.
# Shell commands: Step1. get columns1-3 using "cut -f1-4", e.g. "Zv9_NA101	10001	20000"
# Step2. count the number of duplicates, remove them and add a unique counter (first tab area) using "uniq --count",
# e.g. "	2 Zv9_NA101	10001	20000". The number of duplicates == number of elements in each window
# Step3. Sort the number of duplicates in reverse order (based on the unique counter) using "sort -k1rn"

try:
    procOverlapSelect = subprocess.Popen("bash -c \"overlapSelect -mergeOutput -selectFmt=bed -inFmt=bed "
                                         + outbedSort.name + " " + winCompilebed.name
                                         + " stdout | cut -f1-3 | uniq --count | sort -k1rn > "
                                         + noOfElementCount.name + "\"", shell=True)
    procOverlapSelect.wait() # wait for procOverlapSelect to terminate
except:
    print "Can't open file noOfElementCount"

winCompilebed.close()
outbedSort.close()
noOfElementCount.seek(0)

# check if noOfElementCount has any entries, i.e. there is at least one window containing element(s)
# reading noOfElementCount file and convert each line into a tuple within a master list
if os.stat(noOfElementCount.name).st_size != 0:  #  Check if filesize != 0
    win_list = [] # master list of windows sorted according to noOfElementCount
    for line in noOfElementCount:
        if not line.strip():  # skip blank lines
            continue
        else:
            winLine = line.split()
            winElement = (int(winLine[0]), winLine[1], int(winLine[2]), int(winLine[3]))  
            # with winLine[0] = number of elements, winLine[1] = chrom, winLine[2] = chromStart, winLine[3] = chromEnd
            win_list.append(winElement)

    print "\nNumber of windows that contain elements:", len(win_list)
    if len(win_list) > 1:  # need at least 2 windows containing elements to do statistics.
    # do statistics for number of elements
        sortNumberOfElement_list = [j[0] for j in win_list]  # create list by taking the first item of each tuple in winList
        sortNumberOfElement_list.sort()  # sorts list (from a descending order) to ascending order
        winStatq1 = percentile(sortNumberOfElement_list, 0.25)  # calling percentile function to calculate q1
        winStatmed = percentile(sortNumberOfElement_list, 0.50)
        winStatq3 = percentile(sortNumberOfElement_list, 0.75)
        print "Statistics for number of elements in sliding windows (window size:%s, offset: %s):" % (args.winSize, args.offset)
        winMean = average(sortNumberOfElement_list)  # calling average function
        print "q1: %-12smed: %-12sq3: %-12smean: %-12.3f min: %-12s max: %-12s " % (winStatq1, winStatmed, winStatq3,
                                                                                    winMean, sortNumberOfElement_list[0], sortNumberOfElement_list[-1])
    noOfElementCount.close()

    # show the top 5 unique clusters (most elements overlapping). The tuples in winList is already sorted from the highest number of elements.
    # windows that overlap windows with higher counts are skipped in order to show the top5 non-overlapping windows
    cluster_list = []  # list of selected clusters
    nonOverlappingClusters = isElementOverlappingWithList(win_list, cluster_list) # call function to check if an element is overlapping with other elements in a list
    print "Top 5 non-overlapping clustering windows:(chr\tchrStart\tchrEnd\tnoOfElement)"
    for i, j in zip(range(5), cluster_list):
        print "%s\t%s\t%s\t%s" % (j[1], j[2], j[3], j[0])

    if args.verbose:
        print "-v Top 15 overlapping windows for comparison"
        for i, j in zip(range(15), win_list): 
            print "-v %s\t%s\t%s\t%s" % (j[1], j[2], j[3], j[0])

    # create outfile if user invokes '-winOut [filename]'
    regexwinOutNum = re.compile(r'([><=]+)([0-9]+)')
    winMatch = regexwinOutNum.search(args.winOutNum)
    if args.winOut:
        if winMatch:  # if match exists
            winCondition = winMatch.group(1)  # retrieve conditional operator as winCondition
            winNumber = winMatch.group(2)  # retrieve window size number as winNumber
            winNumber = int(winNumber)  # converts winNumber to integer
            if winNumber >= 1:  #  checks minimum value for winNumber
                print "-winOut file created:", args.winOut.name, "; Parameter -winOutNum:", args.winOutNum
                winNumber = condition(winCondition, winNumber)  # call condition function
                if winCondition == '<' or winCondition == '<=':
                    for i, j in zip(range(len(win_list)), win_list):
                        if j[0] < winNumber:
                            args.winOut.write("%s\t%s\t%s\t%s\n" % (j[1], j[2], j[3], j[0]))
                            verboseprint('selected window: index, winChrom, winStart, winEnd, numberOfElement:',
                                         i, j[1], j[2], j[3], j[0])
                elif winCondition == '>' or winCondition == '>=':
                    for i, j in zip(range(len(win_list)), win_list):
                        if j[0] > winNumber:
                            args.winOut.write("%s\t%s\t%s\t%s\n" % (j[1], j[2], j[3], j[0]))
                            verboseprint('selected window: index, winChrom, winStart, winEnd, numberOfElement:',
                                         i,j[1], j[2], j[3], j[0])
                args.winOut.close()
        else:
            print "Incorrect input value for parameter -overlapOutDist. Output file process aborted."

else:  # if size noOfElementCount.name == 0, ie. no windows were found to contain any elements
    print "No overlaps found between elements and windows. Try reducing window size and offset."

if args.verbose:
    print '-v *finished section on clustering/windows'

# summary of created output files
print "\nSummary of output files created:"
fileCreated = False
if args.spacingOut:
    if os.stat(args.spacingOut.name).st_size != 0:  #  Check if filesize != 0
        print "-spacingOut file created: ", args.spacingOut.name, "; Parameter -spacingOutDist:", args.spacingOutDist
        fileCreated = True
    else:
        print "Empty -spacingOut file or file not created: ", args.spacingOut.name
        print "\tNo spacing elements found within the specified values or incorrect values given for -spacingOutDist:", args.spacingOutDist
        print "\tAccepted values for -spacingOut: \'>=1\'"

if args.overlapOut:
    if os.stat(args.overlapOut.name).st_size != 0:  #  Check if filesize != 0
        print "-overlapOut file created: ", args.overlapOut.name, "; Parameter -overlapOutDist:", args.overlapOutDist
        fileCreated = True
    else:
        print "Empty -overlapOut file or file not created: ", args.overlapOut.name
        print "\tNo overlapping elements found within the specified values or incorrect values given for -overlapOutDist:", args.overlapOutDist
        print "\tAccepted values for absolute overlap: \'a>=1\'; relative overlap: \'0< r<= 1.0\'"

if args.winOut:
    if os.stat(args.winOut.name).st_size != 0:
       print "-winOut file created: ", args.winOut.name, "; Parameter -winOutNum: ", args.winOutNum
       fileCreated = True
    else:
       print "Empty -winOut file or file not created: ", args.winOut.name
       print "\tNo windows found within the specified values or incorrect values given for -winOutNum:", args.winOutNum
       print "\tAccepted values for -winOutNum: \'>=1\'"

if args.keepTmpFiles:
    print "tempFiles saved:"
    print "\toutbedSort logged in:", outbedSort.name
    print "\tbedCov output logged in:", outbedCov.name
    print "\tList of generated windows logged in:", winCompilebed.name
    print "\toverlapSelect output logged in:", noOfElementCount.name
    fileCreated = True

if not fileCreated:
    print "None"

# show histograms, use at end to avoid interrupting script
if args.noHist:  # if true, ie. args.noHist was not called
    plt.show()

################### End script ################################
