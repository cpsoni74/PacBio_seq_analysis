#Test Main 

"""
Edit annotations on the genbank file
Generate an sam alignment file using minimap2
Sort the sam file and generate a bam and bai file so that its all indexed and thereby quicker to parse through
Filter barcodes based on alignment and Q-score & store these reads in a temp file, sorted and indexed
Use the temp file to get corresponding aarS sequences and generate a dataframe
Plot barcodes and mutations
"""

import os
import tempfile
import warnings
import Bio.SeqIO
from Bio import SeqIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
import os, io, random
import pandas as pd
import numpy as np
import pysam
from pysam import bcftools
from thefuzz import fuzz
from thefuzz import process
import matplotlib
import matplotlib.pyplot as plt
#import seaborn as sns
import time
import math
from Read_gb_file import read_template
from gb_feature_edit_parameters import annotating_gb
from PBQscore import qscore
from histo_bc_qscore import histo
from aaRS_qscore import extract_aaRS
from mut_Dict import gen_mut_table
from bcpermut import barcodes_per_mutation
from aaRSperbc import aaRS_per_barcode
""" 
---------------------------------------------------------------------------------
Load in .gb file and use it to find relevant features and locations with location\
---------------------------------------------------------------------------------

---------------------------------------------------------------------------------
#Use genbank file to get the barcode and library site annotation. 
#Open genbank file in txt, note down features and locations to edit as follows
#Modify genbank file annotation with edit_gb_file to replace misc_feature with barcode. Simply open genbank as text and see which feature you need to edit. If it is 2nd, you will use index 1 as counting starts from 0
#example, this is the function definition
#def change_feature_type(template, i, strtloc, endloc, new_type):
#template is gb file, i is the index of the feature, strtloc is the start of the feature as mentioned in gb file, endloc is end, new_type is the new annotation/new feature type
#Open gb_feature_edit_parameters and input any number of changes you want to make
---------------------------------------------------------------------------------

---------------------------------------------------------------------------------
#Input the location for genbank template and also the sam file
#Add bam file and bai file location (only the directory for these files has to exist, the code will generate these files)
---------------------------------------------------------------------------------
"""
Dict = {}

template = '/Users/chintansoni/Desktop/MaPyl Barcoding library seq/PacBio input/Library 1_ L125, N166, V168 NNK/nnk-library-aligned-with-pacbio-reads.gb'
input_sam = '/Users/chintansoni/Desktop/MaPyl Barcoding library seq/PacBio input/Library 1_ L125, N166, V168 NNK/CS_NNK_pb.sam'
input_bam = '/Users/chintansoni/Desktop/MaPyl Barcoding library seq/PacBio input/Library 1_ L125, N166, V168 NNK/CS_NNK_pb_pysamsorted.bam'
bai_file = '/Users/chintansoni/Desktop/MaPyl Barcoding library seq/PacBio input/Library 1_ L125, N166, V168 NNK/CS_NNK_pb_pysamsorted.bam.bai'

change_gb = True #Always set to true 

sort_sam = False
index_file = False
barcode_filter = True #Check qscore in PBQscore
fitlered_Reads_file = True #sort and index filtered barcode reads
aaRS_filter = True #filter_aaRS reads
qcutoff = 20 #qscore cutoff for aaRS reads
mutations = True #generate a dataframe of aaRS mutations and corresponding barcodes
bc_mut = True #generate a dataframe for number of barcodes per unique mutant
mut_bc = True #generate a dataframe for number of unique mutants per barcode

temp_rec = read_template(template)

if change_gb :
    template = annotating_gb(template)   #Opens gb_features_edit_parameters and makes changes to gb and returns edited template
    print ("\n updated file : ")
    print (str(template))
    temp_rec = read_template(template)

print ("\n Features : ")
print (temp_rec[1])

#Run minimap2 after editing features
# for the sam file, use minimap2 to create consensus sequence. Use fasta file for template
#install minimap2 via pip
#minimap2 --MD -Lax map-pb template_file.fasta fastq_data.fastq > file_to_create.sam
#sort the sam file and generate bam file using the following pysam command

if sort_sam:
    pysam.sort("-o", "/Users/chintansoni/Desktop/MaPyl Barcoding library seq/PacBio input/Library 1_ L125, N166, V168 NNK/CS_NNK_pb_pysamsorted.bam","/Users/chintansoni/Desktop/MaPyl Barcoding library seq/PacBio input/Library 1_ L125, N166, V168 NNK/CS_NNK_pb.sam")

#Index the bam file 

if index_file:
    pysam.index(input_bam)

#Filter barcodes based on quality score and alignment and generate a histogram of barcodes
#opens PBQscore.py
#Make sure to enter total number of reads tot and length of barcode bclen in PBQscore.py

if barcode_filter:
    qscore(temp_rec, input_bam, bai_file)
    final_BCDF = pd.read_csv('barcodecount_qscore.csv')
    qual = pd.read_csv('barcode_qscore.csv')
    histo (final_BCDF, qual)

#Collect the filtered reads and generate a bam and bai file

if fitlered_Reads_file:
    pysam.sort("-o", "tempReadssorted.bam", "tempReads.bam")
    pysam.index("tempReadssorted.bam")

#Filter aaRS

if aaRS_filter:

    extract_aaRS(temp_rec, qcutoff)

#Dataframe of aaRS mutations and barcode sequences
#Define MaPylRS AA and MaPylRS DNA sequence in the mut_Dict.py file

if mutations:
    gen_mut_table(template)

#Dataframe & graph of number of barcodes per mutation

if bc_mut:
    barcodes_per_mutation()

#Dataframe and graph of number of unique mutants per barcode

if mut_bc:
    aaRS_per_barcode()

