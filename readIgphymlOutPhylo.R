#!/usr/bin/env Rscript

# expected parameters, in order:
## input igphyml-pass.tab file
## clone group, starting from 1 (largest) to X number of groups available
## output newick file

# E.g usage
# Rscript readIgphymlOutPhylo.R PTFL13_HLP19_1000_igphyml-pass.tab /
# 2 out_PTFL13_2_newick

library(alakazam)
library(ape)

args <- commandArgs(trailingOnly = TRUE)

# read in input file
db = readIgphyml(args[1],format="phylo")

# output selected tree of a clone group in newick format
tree = db$trees[[as.integer(args[2])]]
ape::write.tree(tree, file=args[3])
