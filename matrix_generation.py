import numpy as np
import argparse
import os
from SCOPEseq2_matrix import total_matrix

def parse_user_input():
	parser = argparse.ArgumentParser()
	parser.add_argument('-r','--run-name',required=True,help='Name of sequencing run.')
	parser.add_argument('-data','--data-dir',required=True,help='Name of I/O directory.')
	parser.add_argument('-gene','--genefile',required=True,help='Name of gene reference file.')
	parser.add_argument('-th','--cts_thresh',required=True,help='Molecule threshold.')
	return parser

parser = parse_user_input()
user_input = parser.parse_args()

statsfile = user_input.data_dir+'/'+user_input.run_name+'.exon_address.stats.txt'
genefile = user_input.genefile
filtercts_file = user_input.data_dir+'/'+user_input.run_name+'.exon_address.filt.txt'
matrix_file = user_input.data_dir+'/'+user_input.run_name+'.matrix.txt'
hist_file = user_input.data_dir+'/'+user_input.run_name+'.matrix.hist.txt'
cts_thresh = int(user_input.cts_thresh)

total_matrix(statsfile,genefile,filtercts_file,matrix_file,hist_file,cts_thresh)


