#! /usr/bin/python
import numpy as np

def total_matrix(statsfile,genefile,filtercts_file,matrix_file,hist_file,cts_thresh):
	barcodes = {}
	bclist = []
	i=0
	with open(statsfile) as f:
		for line in f:
			llist = line.split()
			cts = int(llist[2])
			if cts > cts_thresh:
				bc = llist[0]
				bclist.append(bc)
				barcodes[bc] = i
				i+=1
			else:
				break
	gidgene = {}
	matrix = {}
	with open(genefile) as f:
		for line in f:
			llist = line.split()
			gid = llist[0]
			gidgene[gid] = llist[1]
			matrix[gid] = [0 for bc in bclist]	
	with open(filtercts_file) as f:
		for line in f:
			llist = line.split()
			bc = llist[0]
			if bc in barcodes.keys():
				matrix[llist[2]][barcodes[bc]]+=1
	with open(matrix_file,'w') as g:
		st = 'gids\tgenes\t'+'\t'.join(bclist)+'\n'
		g.write(st)
		for gid in gidgene.keys():
			gene = gidgene[gid]
			st = gid+'\t'+gene+'\t'+'\t'.join([str(pt) for pt in matrix[gid]])+'\n'
			g.write(st)
	matrix = np.array(list(matrix.values()))
	norm = sum(matrix)
	with open(hist_file,'w') as g:
		for bc,n in zip(bclist,norm):
			st = bc+'\t'+str(n)+'\n'
			g.write(st)
	return 0
