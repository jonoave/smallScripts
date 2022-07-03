# smallScripts
Personal collection of small tools to manipulate sequence files

## A. General utility
1. *checkmd5sum.py* : Checks md5sum generated from a file vs md5sum provided for a file (missing lines, WIP).
2. *renameFilesRecursive.py*: Based on an input file of old/new names, will rename all files recursively in an input directory. This was designed for file prefixes - to rename random prefixes with QBiC barcodes.   
3. *matchColumnInOneFile.py*: Sort by column1 in file2 based on the order of column1 in file1. Requires pandas package.


## B. Manipulating fasta files (primarily phylogenetic purposes)
1. *format_fasta_oneline.py* : Converts a multi-line sequence fasta file into a (seq1Name, newline, seq1 in one line, newline, seq2Name..). Just uses awk, so no dependencies biopython or perl.
2. *filterMinSpecies.py* : Goes through each fasta file in an input directory, checks the number of sequences in a file, and only copies file to the output directory if number of sequences >= `minNumSpecPer` on `totalSpecies`. E.g if `totalSpecies` is 60 and `minNumSpecPer` is 50, a fasta file must have >30 sequences to be copied to output folder.
3. *seqLengthDist.py*: #  1)get a distribution count of sequence length and 2)generate a histogram of length distribution.
4. *statsFasta.py*: output stats for a fasta file (number of seqs, average seq length, longestSeq, shortestSeq, %GC)

## C. Utility scripts of IgPhyML (B cell phylogenetic inference package)
1. getIgphymlMSAclone.py: Parse the input file for `igPhyML`, typically named xxx_db-pass_productive-T_clone-pass_germ-pass.tsv to create a standard MSA file in fasta format.
2. *readIgphymlOutPhylo.R* : Rscript template to read in output from `igPhyML`, typically xxx_igphyml-pass.tab, to generate trees in newick format.