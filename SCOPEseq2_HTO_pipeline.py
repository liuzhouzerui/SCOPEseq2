#! /usr/bin/python
from SCOPEseq2_demultiplexer import get_cbc_umi, get_csbc
from SCOPEseq2_address import get_cs_address
from SCOPEseq2_addressct import cs_addressct, umifilter
import sys
import numpy as np
import argparse
import os

def parse_user_input():
	parser = argparse.ArgumentParser()
	parser.add_argument('-r','--run-name',required=True,help='Name of sequencing run.')
	parser.add_argument('-data','--data-dir',required=True,help='Name of I/O directory.')
	parser.add_argument('-r1','--read1-fastq',required=True,help='Path to gzipped read 1 fastq file.')
	parser.add_argument('-r2','--read2-fastq',required=True,help='Path to gzipped reqd 2 fastq file.')
	parser.add_argument('-bc1','--barcode1-infile',required=True,help='Path to file containing first cell-identifying barcode segment.')
	parser.add_argument('-bc2','--barcode2-infile',required=True,help='Path to file containing second cell-identifying barcode segment.')
	parser.add_argument('-ind','--index-infile',required=True,help='Path to file containing the Truseq index sequences.')
	parser.add_argument('-b','--citeseq-barcodes',required=True,help='Path to 2-column file with cite-seq features and barcodes.')
	return parser

parser = parse_user_input()
user_input = parser.parse_args()

# Read in cell-identifying barcodes
barcode1_list = [line.split()[0] for line in open(user_input.barcode1_infile)]
barcode2_list = [line.split()[0] for line in open(user_input.barcode2_infile)]
index_list = [line.split()[0] for line in open(user_input.index_infile)]

# Create log file
logfile = user_input.data_dir+'/'+user_input.run_name+'.log'
log = open(logfile,'w')

# Get cell-identifying barcodes, indices, and UMIs from read 1 and error-correc the cell-identifying barcodes and indices
demux_file =  user_input.data_dir+'/'+user_input.run_name+'.citeseq_demux_R1.txt'
demux_N,reads_N = get_cbc_umi(barcode1_list,barcode2_list,index_list,user_input.read1_fastq,demux_file) # get cell-identifying barcodes (CBCs) and unique molecular identifiers (UMIs)
demux_P = float(demux_N)/float(reads_N)*100. # percentage of reads demultiplexed
print('Found %(demux_N)d reads with CBC/UMI out of %(reads_N)d reads or %(demux_P)f%% demultiplexed...' % vars())
log.write('Found %(demux_N)d reads with CBC/UMI out of %(reads_N)d reads or %(demux_P)f%% demultiplexed.\n' % vars())

# Get citeseq feature from read 2 (contains read id, feature information)
fq2file = user_input.read2_fastq
citeseqfile = user_input.citeseq_barcodes
citeseq_demux_file = user_input.data_dir+'/'+user_input.run_name+'.citeseq_demux_R2.txt'
demux_N,reads_N = get_csbc(fq2file,citeseqfile,citeseq_demux_file)
demux_P = float(demux_N)/float(reads_N)*100. # percentage of reads demultiplexed
print('Found %(demux_N)d reads with cite-seq feature out of %(reads_N)d reads or %(demux_P)f%% demultiplexed...' % vars())
log.write('Found %(demux_N)d reads with cite-seq feature out of %(reads_N)d reads or %(demux_P)f%% demultiplexed.\n' % vars())

# Get citeseq feature address from read 2 (contains read id, index, cell-identifying barcode, UMI, feature information)
addressfile = user_input.data_dir+'/'+user_input.run_name+'.citeseq_address.txt.gz'
feature_N,reads_N = get_cs_address(demux_file,citeseq_demux_file,addressfile)
feature_P = float(feature_N)/float(reads_N)*100.
print('Found %(feature_N)d reads with CBC/UMI and cite-seq feature out of %(reads_N)d reads or %(feature_P)f%% cite-seq decoded...' % vars())
log.write('Found %(feature_N)d reads with CBC/UMI and cite-seq feature out of %(reads_N)d reads or %(feature_P)f%% cite-seq decoded.\n' % vars())

# Collapse UMIs
addressct_file = user_input.data_dir+'/'+user_input.run_name+'.citeseq_address.fcts.txt'
molec_N = cs_addressct(addressfile,addressct_file)
print('Found %(molec_N)d cite-seq molecules from %(feature_N)d cite-seq reads...' % vars())
log.write('Found %(molec_N)d cite-seq molecules from %(feature_N)d cite-seq reads...' % vars())  

# Error-correction of UMIs for cite-seq molecules
addressfilt_file = user_input.data_dir+'/'+user_input.run_name+'.citeseq_address.filt.txt'
filt_N, init_N = umifilter(addressct_file,addressfilt_file)
filt_P = (1.0-float(filt_N)/float(init_N))*100.
print('Found %(filt_N)d molecules after filtering %(init_N)d unfiltered molecules for a filter rate of %(filt_P)f%%...' % vars())
log.write('Found %(filt_N)d molecules after filtering %(init_N)d unfiltered molecules for a filter rate of %(filt_P)f%%...' % vars())

log.close()
