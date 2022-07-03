# misc_util
Personal collection of small tools to manipulate sequence files

## A. General utility
1. *checkmd5sum.py*
   
2. *renameFilesRecursive.py*
3. *matchColumnInOneFile.py*
4. 


## B. Manipulating fasta files 
1. *format_fasta_oneline.py* : Converts a multi-line sequence fasta file into a (seq1Name, newline, seq1 in one line, newline, seq2Name..). Just uses awk, so no dependencies biopython or perl.
2. *seqLengthDist.py*
3. *statsFasta.py*

## C. Manipulating sequences in fasta files (primarily phylogenetic purposes)
2. *filterMinSpecies.py* : Goes through each fasta file in an input directory, checks the number of sequences in a file, and only copies file to the output directory if number of sequences >= `minNumSpecPer` on `totalSpecies`. E.g if `totalSpecies` is 60 and `minNumSpecPer` is 50, a fasta file must have >30 sequences to be copied to output folder.
3. *getLongestSeqPerSp.py* : 
4. *seqLengthDist.py*


## D. Utility scripts of IgPhyML (B cell phylogenetic inference package)
1. getIgphymlMSAclone.py
2. readIgphymlOutPhylo.R