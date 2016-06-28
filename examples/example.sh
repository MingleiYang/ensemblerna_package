#!/bin/bash

# This file includes example uses for EnsembleRNA. You can test them by copy and pasting them into
# the command line or by running this script. These are only a few examples to help you get started.
# Try mixing and matching the different commands to personalize your visualization.


#######################################################################################################
# Show the help for EnsembleRNA
#######################################################################################################
ensemblerna -h


#######################################################################################################
# Run with reference fasta file (map will be created from the reference sequence)
#######################################################################################################
ensemblerna ref_WT.fa example_output


#######################################################################################################
# Run with reference fasta file (map will be created from the reference sequence)
# Ex. Guide the structure prediction using SHAPE data
#######################################################################################################
ensemblerna ref_WT.fa example_output -sh ref_WT.shape


#######################################################################################################
# Run with reference fasta file and map fasta file (map will be created from the map sequence)
# Ex. Use the wild-type to create a single map and compare different variants as references
#######################################################################################################
ensemblerna ref_WT.fa example_output
ensemblerna ref_MT.fa example_output -m ref_WT.fa


#######################################################################################################
# Run with reference db file (map will be created from reference sequence)
# Ex. Visualize an ensemble of structures generated using a different sampling algorithm
#######################################################################################################
ensemblerna ref_WT.fa example_output -d ref_WT.db


#######################################################################################################
# Run with reference fasta file and map db file (map will be created from dot-bracket structures)
# Ex. Quickly compare mutants without recreating the map dot-bracket each time
#######################################################################################################
ensemblerna ref_WT.fa example_output
ensemblerna ref_MT.fa example_output -md example_output/ref_WT_map.db


#######################################################################################################
# Run with reference fasta file (map will be created from reference sequence)
# Ex. Increase the structural space explored by the map for longer RNAs
#######################################################################################################
ensemblerna ref_WT.fa example_output/ -s 15




