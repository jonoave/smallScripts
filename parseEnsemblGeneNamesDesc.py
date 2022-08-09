#! /usr/bin/python
# Jun Hoe, Lee (2016)
# parse the downloaded Ensembl biomart file to generate a clean output of unique gene IDs
# input file has: 
# ""Ensembl Gene ID	Associated Gene Name	Phenotype description	Description
# ENSG00000210049	MT-TF	Parkinsonism/MELAS overlap syndrome	mitochondrially encoded tRNA phenylalanine [Source:HGNC 
# Symbol;Acc:HGNC:7481]
# ENSG00000210049	MT-TF	MERRF	mitochondrially encoded tRNA phenylalanine [Source:HGNC Symbol;Acc:HGNC:7481]
# ENSG00000211459	MT-RNR1	MERRF	mitochondrially encoded 12S RNA [Source:HGNC Symbol;Acc:HGNC:7470]""
# output should combine genes with multiple phenotype entries into one

import sys

#################################
# A. read input parameters


inputFile = open(sys.argv[1], "r")
outputFile = open(sys.argv[2], "w")    

#################################
# B. parse input file

# create dictionaries
geneIDName_dict = {}
geneIDPhenotype_dict = {}
geneIDdesc_dict = {}


readInput = inputFile.readlines()
for line in readInput[1:]: # skip header line
    try:
        print "line", line
        splitLine = line.split("\t")
        geneID = splitLine[0].replace("\n", "")
        geneName = splitLine[1].replace("\n", "")
        if not geneID in geneIDName_dict:
            geneIDName_dict[geneID] = geneName
        if splitLine[3]: # gene description exists
            geneDesc = splitLine[3].replace("\n", "")
            print "geneDesc", geneDesc
            if not geneID in geneIDdesc_dict:
                geneIDdesc_dict[geneID] = geneDesc
        if splitLine[2]: # gene phenotype exists
            genePhenotype = splitLine[2].replace("\n", "")
            print "gene Phenotype", genePhenotype
            if not geneID in geneIDPhenotype_dict:
                 geneIDPhenotype_dict[geneID] = [genePhenotype]
                 print "genePhenotype entry created in geneIDPhenotype_dict"
            else:
                 geneIDPhenotype_dict[geneID].append(genePhenotype)
                 print "genePhenotype entry appended in geneIDPhenotype_dict"
    except:
        continue

print "geneIDName_dict", geneIDName_dict.items()[:5]
print "geneIDPhenotype_dict",  geneIDPhenotype_dict.items()[:5]
print "geneIDdesc_dict", geneIDdesc_dict.items()[:5]

inputFile.close()

###################################
# C. write out dictionary contents to output file

for key, value in geneIDName_dict.items():
     description = geneIDdesc_dict.get(key, "None") # returs "None" if key doesn't exist
     if not description: # if entry exists, but as an empty string
         description = "None"
     phenotype_list = geneIDPhenotype_dict.get(key) # returns None if key doesn't exist
     if phenotype_list:
         splitPhenotype = str(phenotype_list).strip('[]')
     else: 
         splitPhenotype = "None"
     writeLine = (key + "\t" + value + "\t" + description + "\t"+ splitPhenotype + "\n")
     outputFile.write(writeLine)

outputFile.close()


