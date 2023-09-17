#! /usr/bin/python
import sys
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from collections import Counter
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--filter-file',required=True,help='Name of file containing filtered molecules (e.g. *exon_address.cfilt.txt)')
    parser.add_argument('-r','--reference-file',required=True,help='Two-column file with gids and gene symbols')
    parser.add_argument('-p','--prefix',required=True,help='Prefix for output files including path')
    parser.add_argument('-l','--lower-thresh',required=True,help='Lower count threshold parameter for emptyDrops')
    return parser

parser = parse_user_input()
user_input = parser.parse_args()

cfilt_INFILE = user_input.filter_file 
ref_INFILE = user_input.reference_file
output_PREFIX = user_input.prefix
lower_THRESH = int(user_input.lower_thresh)
cfilt_OUTFILE = cfilt_INFILE+'.conv.txt'

cbcs = {}
genes = {}
i=1
j=1
with open(cfilt_OUTFILE,'w') as h:
	with open(cfilt_INFILE) as f:
		for line in f:
			llist = line.split()
			cbc = llist[0]
			gene = llist[2]
			if cbc in cbcs.keys():
				c = cbcs[cbc]
			else:
				c = i
				cbcs[cbc] = i
				i+=1
			if gene in genes.keys():
				g = genes[gene]
			else:
				g = j
				genes[gene] = j
				j+=1
			h.write('%(c)d\t%(g)d\n' % vars())

Rcmd1 = 'library("DropletUtils");'
Rcmd2 = 'data<-read.csv("%(cfilt_OUTFILE)s",header=FALSE,sep="\t",stringsAsFactors=FALSE);' % vars()
Rcmd3 = 'mat<-makeCountMatrix(data[,2],data[,1]);'
Rcmd4 = 'edrops<-emptyDrops(mat,lower=%(lower_THRESH)d,test.ambient=TRUE);' % vars()
Rcmd5 = 'write.table(edrops,file="%(output_PREFIX)s.edrops.txt",sep="\t",col.names=TRUE,row.names=TRUE);' % vars()

ro.r(Rcmd1)
ro.r(Rcmd2)
ro.r(Rcmd3)
ro.r(Rcmd4)
ro.r(Rcmd5)

cbcs_inv = {i:cbc for cbc,i in cbcs.items()}
edrops_OUTFILE = output_PREFIX+'.edrops.txt'
ambient_pvals = []
keep_cbcs = {}
keep = []
totals = []
with open(edrops_OUTFILE) as f:
    next(f)
    for line in f:
        llist = line.split()
        total = int(llist[1])
        if total > lower_THRESH:
            if llist[4] == 'TRUE':
                totals.append(total)
                keep.append(1)
                keep_cbcs[cbcs_inv[int(llist[0][1:-1])]] = []
            else:
                totals.append(total)
                keep.append(0)
        else:
            totals.append(total)
            keep.append(0)
            ambient_pvals.append(float(llist[3]))


with open(cfilt_INFILE) as f:
    for line in f:
        llist = line.split()
        cbc = llist[0]
        if cbc in keep_cbcs.keys():
            keep_cbcs[cbc].append(llist[2])

gdict = {}
with open(ref_INFILE) as f:
    for line in f:
        llist = line.split()
        gdict[llist[0]] = llist[1]

matrix = []
gcts = []
mcts = []
for cbc in keep_cbcs.keys():
    cts = Counter(keep_cbcs[cbc])
    matrix.append([cts[gid] if gid in cts.keys() else 0 for gid in gdict.keys()])
    gcts.append(len(set(keep_cbcs[cbc])))
    mcts.append(len(keep_cbcs[cbc]))
matrix = np.transpose(np.array(matrix))

matrix_OUTFILE = output_PREFIX+'.edrops.matrix.txt'
with open(matrix_OUTFILE,'w') as g:
    for gid,vec in zip(gdict.keys(),matrix):
        gene = gdict[gid]
        st = gid+'\t'+gene+'\t'+'\t'.join([str(pt) for pt in vec])+'\n'
        g.write(st)

hist_OUTFILE = output_PREFIX+'.edrops.matrix.hist.txt'
with open(hist_OUTFILE,'w') as g:
    for cbc,mct,gct in zip(keep_cbcs.keys(),mcts,gcts):
        st = cbc+'\t'+str(mct)+'\t'+str(gct)+'\n'
        g.write(st)

pdf_OUTFILE = output_PREFIX+'.edrops.pdf'
with PdfPages(pdf_OUTFILE) as pdf:
    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(1,2,1)
    total = np.array(totals)
    keep = np.array(keep)
    ind = np.argsort(total)
    total = total[ind][::-1]
    keep = keep[ind][::-1]
    ax1.scatter(range(1,len(total)+1),total,c=keep,cmap='bwr',s=3)
    ax1.set_xlabel('Cell Barcodes')
    ax1.set_ylabel('Number of Molecules Detected')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2 = fig.add_subplot(1,2,2)
    chist,cedge = np.histogram(np.array(ambient_pvals),bins=100)
    ax2.bar(cedge[0:len(cedge)-1],chist,width=0.01,color='k')
    ax2.set_xlabel('p-value')
    ax2.set_ylabel('Number of Ambient Cell Barcodes')
    pdf.savefig()


        


