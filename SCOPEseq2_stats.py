#! /usr/bin/python
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def chist(arr):
	ind = np.argsort(arr)
	srt_arr = arr[ind][::-1]
	tot = float(sum(arr))
	ch = []
	c = 0.0
	for pt in srt_arr:
		c+=float(pt)/tot
		ch.append(c)
	return np.array(ch)

def get_stats(filtercts_file,stats_file):
	cell_genes = {}
	cell_umi = {}
	cell_reads = {}
	with open(filtercts_file) as f:
		for line in f:
			llist = line.split()
			cbc = llist[0]
			if cbc in cell_umi.keys():
				cell_umi[cbc].append(llist[1])
				cell_genes[cbc].append(llist[2])
				cell_reads[cbc].append(int(llist[3]))
			else:
				cell_umi[cbc] = [llist[1]]
				cell_genes[cbc] = [llist[2]]
				cell_reads[cbc] = [int(llist[3])]
	cbclist = list(cell_genes.keys())
	Ngenes = np.array([len(set(cell_genes[cbc])) for cbc in cbclist])
	Numis = np.array([len(cell_umi[cbc]) for cbc in cbclist])
	Nreads = np.array([sum(cell_reads[cbc]) for cbc in cbclist])
	cbclist = np.array(cbclist)
	ind = np.argsort(Numis)
	cbclist = cbclist[ind][::-1]
	Ngenes = Ngenes[ind][::-1]
	Numis = Numis[ind][::-1]
	Nreads = Nreads[ind][::-1]
	with open(stats_file,'w') as g:
		for cbc,ngenes,numis,nreads in zip(cbclist,Ngenes,Numis,Nreads):
			st = cbc+'\t'+str(nreads)+'\t'+str(numis)+'\t'+str(ngenes)+'\n'
			g.write(st)
	return Nreads,Numis,Ngenes

def plot_stats(filtercts_file,stats_file,pdf_file):
	Nreads,Numis,Ngenes = get_stats(filtercts_file,stats_file)
	chist_umis = chist(Numis)
	chist_reads = chist(Nreads)
	rpm = Nreads/Numis
	mpg = Numis/Ngenes

	with PdfPages(pdf_file) as pdf:
		fig = plt.figure(figsize=(10,15))
		ax = fig.add_subplot(3,2,1)
		ax.semilogx(range(1,len(chist_umis)+1),chist_umis,'k')
		ax.set_xlabel('Cell Barcode')
		ax.set_ylabel('Cumulative Histogram for Molecules')
		ax = fig.add_subplot(3,2,2)
		ax.semilogx(range(1,len(chist_reads)+1),chist_reads,'k')
		ax.set_xlabel('Cell Barcode')
		ax.set_ylabel('Cumulative Histogram for Reads')
		ax = fig.add_subplot(3,2,3)
		ax.plot(range(1,len(chist_umis)+1),Numis,'k')
		ax.set_xlabel('Cell Barcode')
		ax.set_ylabel('Number of Molecules')
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax = fig.add_subplot(3,2,4)
		ax.plot(range(1,len(chist_umis)+1),Ngenes,'k')
		ax.set_xlabel('Cell Barcode')
		ax.set_ylabel('Number of Genes')
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax = fig.add_subplot(3,2,5)
		ax.semilogx(range(1,len(chist_reads)+1),rpm,'k')
		ax.set_xlabel('Cell Barcode')
		ax.set_ylabel('Number of Reads per Molecule')
		ax = fig.add_subplot(3,2,6)
		ax.semilogx(range(1,len(chist_reads)+1),mpg,'k')
		ax.set_xlabel('Cell Barcode')
		ax.set_ylabel('Number of Molecules per Gene')
		pdf.savefig()
	return 0
	
	
			
					
			
			
