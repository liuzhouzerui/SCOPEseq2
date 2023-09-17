#! /usr/bin/python
import gzip
import io
from sys import getsizeof
import numpy as np

def get_size(obj, seen=None):
	"""Recursively finds size of objects"""
	size = getsizeof(obj)
	if seen is None:
		seen = set()
	obj_id = id(obj)
	if obj_id in seen:
		return 0
	# Important mark as seen *before* entering recursion to gracefully handle
	# self-referential objects
	seen.add(obj_id)
	if isinstance(obj, dict):
		size += sum([get_size(v, seen) for v in obj.values()])
		size += sum([get_size(k, seen) for k in obj.keys()])
	elif hasattr(obj, '__dict__'):
		size += get_size(obj.__dict__, seen)
	elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
		size += sum([get_size(i, seen) for i in obj])
	return size

# Enumerate HD=1 sequences for a given barcode sequence
def enumerate_bc(bc):
	bcs = []
	N = len(bc)
	for i in range(N):
		bc1=''
		bc2=''
		bc3=''
		bc4=''
		bc5=''
		for j in range(N):
			if i==j:
				bc1+='A'
				bc2+='G'
				bc3+='C'
				bc4+='T'
				bc5+='N'
			else:
				bc1+=bc[j]
				bc2+=bc[j]
				bc3+=bc[j]
				bc4+=bc[j]
				bc5+=bc[j]
		bcs.append(bc1)
		bcs.append(bc2)
		bcs.append(bc3)
		bcs.append(bc4)	
		bcs.append(bc5)
	bcs = list(set(bcs))
	return bcs

# Generate dictionary linking HD=1 sequences (keys) to correct barcode (values)
def get_cbc_dict(bcs):
	cbc_dict = {}
	for bc,i in zip(bcs,range(len(bcs))):
		allbcs = enumerate_bc(bc)
		for allbc in allbcs:
			cbc_dict[allbc] = i
	return cbc_dict	

# Extract SCOPEseq2 cell barcodes (CBCs) and unique molecular identifiers (UMIs) for read 1 fastq
def get_cbc_umi(first_bc,second_bc,index_bc,read1_fastq,demux_file):
	first_bc_dict = get_cbc_dict(first_bc)
	second_bc_dict = get_cbc_dict(second_bc)
	index_bc_dict = get_cbc_dict(index_bc)
	productive=0
	total=0
	with open(demux_file,'w') as g:
		with io.BufferedReader(gzip.open(read1_fastq,'rb')) as f:
			i=0
			for line in f:
				if i == 0:
					dlist = line.decode().split()
					readid = ':'.join(dlist[0].split(':')[3:7])
					index = dlist[1].split(':')[3]
					if index in index_bc_dict.keys():
						index = index_bc_dict[index]
					else:
						index = '-1'
				elif i == 1:
					if index != '-1':
						dline = line.decode()
						bc1 = dline[2:10] # first 8-base barcode segment
						if bc1 in first_bc_dict.keys():
							bc1 = first_bc_dict[bc1]
							bc2 = dline[12:20] # second 8-base barcode segment
							if bc2 in second_bc_dict.keys():
								bc2 = second_bc_dict[bc2] 
								umi = dline[0:2]+dline[10:12]+dline[20:24] # 8-base UMI in 2x2-base and 1x4-base blocks
								if umi.find('N') == -1:
									barcode = str(index)+'_'+str(bc1)+'_'+str(bc2)+'_'+umi
									productive+=1
								else:	
									barcode='-1'
							else:
								barcode='-1'
						else:
							barcode = '-1'
					else:
						barcode = '-1'
				i+=1
				if i==4:
					i=0
					st = readid+'\t'+barcode+'\n'
					g.write(st)
					total+=1
	return productive,total


def get_csbc(fq2file,citeseqfile,citeseq_demux_file):
	# get enumerated barcodes for cite-seq antibodies
	csbc_dict = {}
	with open(citeseqfile) as f:
		for line in f:
			llist = line.split()
			bc = llist[1]
			bcs = enumerate_bc(bc)
			for bc2 in bcs:
				csbc_dict[bc2] = llist[0]
	# decode cite-seq features from read 2 fastq
	cslen = len(list(csbc_dict.keys())[0])
	productive=0
	total=0
	with open(citeseq_demux_file,'w') as g:
		with io.BufferedReader(gzip.open(fq2file,'rb')) as f:
			i=0
			for line in f:
				if i == 0:
					dlist = line.decode().split()
					readid = ':'.join(dlist[0].split(':')[3:7])
				elif i == 1:
					dline = line.decode()
					csbc = dline[0:cslen]  # first 15-base cite-seq barcode segment
					if csbc in csbc_dict.keys():
						feature = csbc_dict[csbc]
						productive+=1
					else:
						feature = '-1'
				i+=1
				if i==4:
					i=0
					st = readid+'\t'+feature+'\n'
					g.write(st)
					total+=1
	return productive,total



	
