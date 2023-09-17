#! /usr/bin/python
import gzip
import io
import numpy as np

def addressct(addressfile,gene_ctfile,exon_ctfile):
	address_dict = {}
	with io.BufferedReader(gzip.open(addressfile,'rb')) as f:
		for line in f:
			dllist = line.decode().split()
			address = "\t".join(dllist[1:4])
			case = int(dllist[4])
			if case == 0: # ambiguous gene, unambiguous exon
				if address in address_dict.keys():
					address_dict[address][1]+=1
				else:
					address_dict[address] = [0,1]
			elif case == 1: # unambiguous gene, unambiguous exon
				if address in address_dict.keys():
					address_dict[address][0]+=1
					address_dict[address][1]+=1
				else:
					address_dict[address] = [1,1]
			elif case == 2: # unambiguous gene, no exon
				if address in address_dict.keys():
					address_dict[address][0]+=1
				else:
					address_dict[address] = [1,0]
	gct = 0
	ect = 0
	with open(gene_ctfile,'w') as g1:
		with open(exon_ctfile,'w') as g2:
			for address in address_dict.keys():
				pt1 = address_dict[address][0]
				pt2 = address_dict[address][1]
				if pt1 > 0:
					st = address+'\t'+str(pt1)+'\n'
					g1.write(st)
					gct+=1
				if pt2 > 0:
					st = address+'\t'+str(pt2)+'\n'
					g2.write(st)
					ect+=1
	return gct,ect

# The purpose of this function is to apply an H=1 Hamming filter to the UMIs
# in a DropSeq experiment. Note that this function is partially converted into
# C code using cython and compiled for speed.
def umifilter(addresscts_file,filtercts_file): # input file is the output of DropSeqPipeline6_addresscts.py
	gene_dict = {} # dictionary where keys are gene symbols, values are tuples of cell barcode, UMI sequence, number of reads
	init_cts = 0
	filt_cts = 0
	with open(addresscts_file) as f:
		for line in f:
			llist = line.split()
			gene = llist[2] # gene symbol
			tup = llist[0],llist[1],int(llist[3]) # tuple containing cell barcode, UMI sequence, and number of reads
			if gene in gene_dict.keys():
				gene_dict[gene].append(tup)
			else:
				gene_dict[gene] = [tup]
			init_cts+=1
	with open(filtercts_file,'w') as g:
		for gene in gene_dict.keys(): # only need to compare UMIs with the same gene symbol
			tups = gene_dict[gene]
			cells = {} # pre-sort based on cell barcodes- make a dictionary where the keys are cell barcodes, values are tuple with UMI sequence, number of reads
			for tup in tups:
				if tup[0] in cells.keys():
					cells[tup[0]].append((np.array(list(tup[1])),tup[2]))
				else:
					cells[tup[0]] = [(np.array(list(tup[1])),tup[2])]
			keep = []
			for cell in cells.keys(): # only need to compare UMIs with the same cell barcode
				tups = cells[cell] 
				for tup in tups:
					umi = tup[0]
					cts = tup[1]
					go = 1
					for tup2 in tups: # loop through all UMI-UMI pairs with the same cell barcode and gene symbol
						if tup2[1] > cts: # if you find a second UMI with more reads than the first UMI
							d = np.count_nonzero(umi!=tup2[0]) 
							if d < 2: # and that second UMI is with a Hamming distance of 1 from the first
								go = 0 # the first UMI is discarded
								break
					if go == 1:
						st = cell+'\t'+''.join(umi)+'\t'+gene+'\t'+str(cts)+'\n'
						g.write(st)
						filt_cts+=1
	return filt_cts,init_cts


def cs_addressct(addressfile,addressct_file):
	address_dict = {}
	cts = 0
	with io.BufferedReader(gzip.open(addressfile,'rb')) as f:
		for line in f:
			dllist = line.decode().split()
			address = "\t".join(dllist[1:4])
			if address in address_dict.keys():
				address_dict[address]+=1
			else:
				address_dict[address]=1
	with open(addressct_file,'w') as f:
		for address in address_dict.keys():
			st = address+'\t'+str(address_dict[address])+'\n'
			f.write(st)
			cts+=1
	return cts
