#! /usr/bin/python
from SCOPEseq2_demultiplexer import get_cbc_umi
from SCOPEseq2_clipper import clipper
from SCOPEseq2_address import get_address
from SCOPEseq2_addressct import addressct, umifilter
from SCOPEseq2_stats import get_stats, plot_stats
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
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
	parser.add_argument('-ind','--index-infile',required=True,help='Path to file containing the index sequences.')
	parser.add_argument('-ref','--reference-dir',required=True,help='Name of reference genome directory for STAR.')
	parser.add_argument('-gtf','--gtf',required=True,help='Path to transcriptome annotation gtf file for STAR.')
	parser.add_argument('-dist','--overhang_distance',default=65,help='Sets the sjdbOverhang parameter for STAR.')
	parser.add_argument('-t','--threads',default=1,help='Sets the number of threads for STAR and samtools.')
	return parser

parser = parse_user_input()
user_input = parser.parse_args()

# Read in cell-identifying barcodes and Illumina indices
barcode1_list = [line.split()[0] for line in open(user_input.barcode1_infile)]
barcode2_list = [line.split()[0] for line in open(user_input.barcode2_infile)]
index_list = [line.split()[0] for line in open(user_input.index_infile)]

# Create log file
logfile = user_input.data_dir+'/'+user_input.run_name+'.log'
log = open(logfile,'w')

# Get cell-identifying barcodes, indices, and UMIs from read 1 and error-correc the cell-identifying barcodes and indices
demux_file =  user_input.data_dir+'/'+user_input.run_name+'.demux.txt'
demux_N,reads_N = get_cbc_umi(barcode1_list,barcode2_list,index_list,user_input.read1_fastq,demux_file) # get cell-identifying barcodes (CBCs) and unique molecular identifiers (UMIs)
demux_P = float(demux_N)/float(reads_N)*100. # percentage of reads demultiplexed
print('Found %(demux_N)d reads with CBC/UMI out of %(reads_N)d reads or %(demux_P)f%% demultiplexed...' % vars())
log.write('Found %(demux_N)d reads with CBC/UMI out of %(reads_N)d reads or %(demux_P)f%% demultiplexed.\n' % vars())

ref = user_input.reference_dir
gtf = user_input.gtf
t = user_input.threads
fq2file = user_input.read2_fastq
fq2clip = user_input.data_dir+'/'+user_input.run_name+'_R2.clip.fastq.gz'
bamfile = user_input.data_dir+'/'+user_input.run_name+'.'
dist = user_input.overhang_distance

# Clip poly(A) tails from read 2 and add barcode information to a new read 2 fastq
clipped_N,total_N = clipper(fq2file,fq2clip)
clipped_P = float(clipped_N)/float(total_N)*100.
print('Found %(clipped_N)d reads with poly(A) out of %(total_N)d reads or %(clipped_P)f%% clipped...' % vars())
log.write('Found %(clipped_N)d reads with poly(A) out of %(total_N)d reads or %(clipped_P)f%% clipped.\n' % vars())

# Two-pass, splice-/annotatin-aware alignment wtih STAR
cmd = 'STAR --readFilesCommand gzip -d -c --genomeDir %(ref)s --sjdbOverhang %(dist)s --sjdbGTFfile %(gtf)s --twopassMode Basic --runThreadN %(t)s --readFilesIn %(fq2clip)s --outFileNamePrefix %(bamfile)s --outSAMtype BAM Unsorted' % vars()
print('STAR command...')
print(cmd)
os.system(cmd)

# Create address file of aligned reads from bam file (contains read id, index, cell-identifying barcode, UMI, gene, and alignment case information 
bamfile = bamfile+'Aligned.out.bam'
addressfile = user_input.data_dir+'/'+user_input.run_name+'.address.txt.gz'
uniqalign_N,bam_N = get_address(gtf,bamfile,addressfile,demux_file)
uniqalign_P = float(uniqalign_N)/float(bam_N)*100.
print('Found %(uniqalign_N)d unique gene alignments out of %(bam_N)d aligned reads or %(uniqalign_P)f%% uniquely aligned...' % vars())
log.write('Found %(uniqalign_N)d unique gene alignments out of %(bam_N)d aligned reads or %(uniqalign_P)f%% uniquely aligned.\n' % vars())

# Collapse UMIs
exon_addressct_file = user_input.data_dir+'/'+user_input.run_name+'.exon_address.cts.txt'
gene_addressct_file = user_input.data_dir+'/'+user_input.run_name+'.gene_address.cts.txt'
gene_molec_N,exon_molec_N = addressct(addressfile,gene_addressct_file,exon_addressct_file)
print('Found %(gene_molec_N)d whole-gene molecules and %(exon_molec_N)d exonic molecules from %(uniqalign_N)d unique gene alignments...' % vars())
log.write('Found %(gene_molec_N)d whole-gene molecules and %(exon_molec_N)d exonic molecules from %(uniqalign_N)d unique gene alignments.\n' % vars())

# Error-correction of UMIs for exon-only aligning molecules
exon_addressfilt_file = user_input.data_dir+'/'+user_input.run_name+'.exon_address.filt.txt'
exon_filt_N,exon_init_N = umifilter(exon_addressct_file,exon_addressfilt_file)
exon_filt_P = (1.0-float(exon_filt_N)/float(exon_init_N))*100.
print('Found %(exon_filt_N)d exonic molecules after filtering %(exon_init_N)d unfiltered molecules for a filter rate of %(exon_filt_P)f%%...' % vars())
log.write('Found %(exon_filt_N)d exonic molecules after filtering %(exon_init_N)d unfiltered molecules for a filter rate of %(exon_filt_P)f%%.\n' % vars())

# Error-correction of UMIs for whole gene body aligning molecules
gene_addressfilt_file = user_input.data_dir+'/'+user_input.run_name+'.gene_address.filt.txt'
gene_filt_N,gene_init_N = umifilter(gene_addressct_file,gene_addressfilt_file)
gene_filt_P = (1.0-float(gene_filt_N)/float(gene_init_N))*100.
print('Found %(gene_filt_N)d whole-gene molecules after filtering %(gene_init_N)d unfiltered molecules for a filter rate of %(gene_filt_P)f%%...' % vars())
log.write('Found %(gene_filt_N)d whole-gene molecules after filtering %(gene_init_N)d unfiltered molecules for a filter rate of %(gene_filt_P)f%%\n.' % vars())

stats_file = user_input.data_dir+'/'+user_input.run_name+'.exon_address.stats.txt'
pdf_file = user_input.data_dir+'/'+user_input.run_name+'.exon_address.stats.pdf'
plot_stats(exon_addressfilt_file, stats_file, pdf_file)

log.close()
