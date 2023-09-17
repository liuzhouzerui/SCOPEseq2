#! /usr/bin/python
import gc
from collections import defaultdict
from bx.intervals.intersection import IntervalTree, Interval
from pysam import AlignmentFile
import gzip

def get_interval_tree_for_chromosome(gtf, chrm):
#   Returns interval trees of genes on the given chromosome for each strand,
#   arranged in a {strand: tree} dictionary
	tree_dict = {strand: IntervalTree() for strand in gtf.keys()}
	for strand, chrm_dict in gtf.items():
		for iv in chrm_dict[chrm]:
	    # assumes iv of form (left, right, gene_id) as specified in get_gtf
			tree_dict[strand].insert_interval(Interval(*iv))	
	return tree_dict


def get_address(gtfile,bamfile,addressfile,demux_file):
	# For both the + and - strands, there is a default dictionary in the gtf metadictionary with chromosome names as keys and values
	# that are lists of tuples (gene_id, left coordinate, right coordinate) indicate gene regions on that chromosome and strand.
	genegtf = {'+':defaultdict(list),'-':defaultdict(list)} # whole-gene metadictionary keyed by '+' and '-' (indicating strand)
	exongtf = {'+':defaultdict(list),'-':defaultdict(list)} # exon-only metadictionary keyed by '+' and '-' (indicating strand)
	chromosomes = []
	with open(gtfile) as g:
		for line in g:
			if line[0] != '#':
				llist = line.split()
				ch = llist[0]   # chromosome
				x1 = int(llist[3]) # left position
				x2 = int(llist[4]) # right position
				d = llist[6]       # strand (+ or -)
				gid = llist[llist.index('gene_id')+1][1:-2] # gene_id
				region = llist[2] # region type
				if region == 'gene':
					genegtf[d][ch].append((x1,x2+1,gid)) # load gene region into whole-gene metadictionary
				elif region == 'exon':
					exongtf[d][ch].append((x1,x2+1,gid)) # load exon region into exon-only metadictionary
				if ch not in chromosomes:
					chromosomes.append(ch)
	gene_trees = {}
	exon_trees = {}
	for ch in chromosomes:
		gene_trees[ch] = get_interval_tree_for_chromosome(genegtf, ch)
		exon_trees[ch] = get_interval_tree_for_chromosome(exongtf, ch)

	store_mm_gene = {}
	store_mm_exon = {}

	barcode_dict = {}
	with open(demux_file) as f:
		for line in f:
			llist = line.split()
			bc = llist[1]
			if bc != '-1':
				barcode_dict[llist[0]] = bc	

	i = 0 
	j = 0
	cts = 0
	old_readid = ''
	with gzip.open(addressfile, 'wb') as g:
		with AlignmentFile(bamfile,'rb') as f:
			for read in f:
				readid = ':'.join(read.qname.split(':')[3:7])
				if readid in barcode_dict.keys():
					bc = barcode_dict[readid].split('_')
					cellbc = '_'.join(bc[0:3])
					umi = bc[3]
					ch = f.getrname(read.reference_id)
					if ch != None and ch in gene_trees.keys():
						flag = read.flag
						p1 = read.reference_start+1
						p2 = read.reference_end+1
						hits = read.get_tag('NH') 
						if flag == 0 or flag == 256:
							strand = '+'
						elif flag == 16 or flag == 272:
							strand = '-'
						gene_interval = gene_trees[ch][strand].find(p1,p2) # try to get gene from + or - whole-gene tree
						exon_interval = exon_trees[ch][strand].find(p1,p2) # try to get gene from + or - exon-only tree
						exon_gids = [iv.value for iv in exon_interval] # how many genes are associated with the exon intervals?
						n_exon_gids = len(set(exon_gids)) # need to count because multiple exons per gene
						n_gene_intervals = len(gene_interval) # only one "gene per gene", don't need to count gene ids for whole-gene
						if n_gene_intervals > 1:
							if n_exon_gids == 1:	# case 0: whole-gene is ambiguous but exon-only is unambiguous
								gid = exon_gids[0]
								if hits == 1:
									cts+=1
									st = readid+'\t'+cellbc+'\t'+umi+'\t'+gid+'\t0\n'
									g.write(st.encode())
								else: # handle multi-mappers- will keep multimappers if there's only one that maps strand-specifically and within a gene
									if readid not in store_mm_exon: 
										store_mm_exon[readid] = gid,cellbc,umi,0
									else:
										store_mm_exon[readid] = -1 # if readid is already in store_mm, then there must be more than one strand-specific alignment, so discard
						elif n_gene_intervals == 1:     
							if n_exon_gids == 1:	# case 1: whole-gene is unambiguous, therefore exon-only is unambiguous if it exists
								gid = exon_gids[0]
								if hits == 1:
									cts+=1
									st = readid+'\t'+cellbc+'\t'+umi+'\t'+gid+'\t1\n'
									g.write(st.encode())
								else:
									if readid not in store_mm_exon:
										store_mm_exon[readid] = gid,cellbc,umi,1
									else:
										store_mm_exon[readid] = -1
							elif n_exon_gids == 0:	# case 2: whole gene is unambiguous and exon-only does not exist
								gid = gene_interval[0].value
								if hits == 1:
									cts+=1
									st = readid+'\t'+cellbc+'\t'+umi+'\t'+gid+'\t2\n'
									g.write(st.encode())
								else:
									if readid not in store_mm_gene:
										store_mm_gene[readid] = gid,cellbc,umi,2
									else:
										store_mm_gene[readid] = -1
				i+=1
		for readid in store_mm_gene.keys(): # write all multi-mappers that are actually unique
			if store_mm_gene[readid] != -1:
				gid = store_mm_gene[readid][0]
				cellbc = store_mm_gene[readid][1]
				umi = store_mm_gene[readid][2]
				case = store_mm_gene[readid][3]
				cts+=1
				st = readid+'\t'+cellbc+'\t'+umi+'\t'+gid+'\t'+str(case)+'\n'
				g.write(st.encode())
		for readid in store_mm_exon.keys(): # write all multi-mappers that are actually unique
			if store_mm_exon[readid] != -1:
				gid = store_mm_exon[readid][0]
				cellbc = store_mm_exon[readid][1]
				umi = store_mm_exon[readid][2]
				case = store_mm_exon[readid][3]
				cts+=1
				st = readid+'\t'+cellbc+'\t'+umi+'\t'+gid+'\t'+str(case)+'\n'
				g.write(st.encode())
	return cts,i
		

def get_cs_address(demux_file,citeseq_demux_file,addressfile):
	cts=0
	i=0
	with gzip.open(addressfile, 'wb') as g:
		with open(demux_file) as f1, open(citeseq_demux_file) as f2:
			for line1,line2 in zip(f1, f2):
				# get cellbc and umi from read1 demux file
				llist1 = line1.split()
				readid = llist1[0]
				if llist1[1] != '-1':
					bc = llist1[1].split('_')
					cellbc = '_'.join(bc[0:3])
					umi = bc[3]
					# get feature from read2 demux file
					llist2 = line2.split()
					readid2 = llist2[0]
					if readid2 != readid:
						print('Error: readid does not match')
						exit()
					if llist2[1] != '-1':
						feature = llist2[1]
						st = readid+'\t'+cellbc+'\t'+umi+'\t'+feature+'\t'+'0'+'\n'
						g.write(st.encode())
						cts+=1
				i+=1
	return cts,i

#	
