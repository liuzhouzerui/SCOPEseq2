#! /usr/bin/python
import pandas as pd
import anndata as anndata
import scanpy as sc
import argparse
import os

def parse_user_input():
	parser = argparse.ArgumentParser()
	parser.add_argument('-data','--data-dir',required=True,help='Directory of input files.')
	parser.add_argument('-p','--prefix',required=True,help='Prefix of input files.')
	parser.add_argument('-n','--hto-number',required=True,help='Number of HTOs used in the experiment')
	parser.add_argument('-ind','--index-map',required=True,help='Path to the file containing index mapping table between HTO and RNA-seq.')
	return parser

parser = parse_user_input()
user_input = parser.parse_args()

data = user_input.data_dir
prefix = user_input.prefix
HTO_num = int(user_input.hto_number)
index_map_file = user_input.index_map

# Read in HTO counts
fcts = pd.read_csv(data+prefix+'_HTO.citeseq_address.filt.txt', sep='\t', header=None)
fcts.columns = ['cell_barcode', 'UMI', 'hashtag', 'count']
fcts_matrix = fcts[['cell_barcode', 'hashtag', 'count']][fcts['count']>=0].groupby(['cell_barcode', 'hashtag']).count().unstack(level=-1, fill_value=0)['count']
fcts_matrix = fcts_matrix.T[fcts_matrix.sum().rank(ascending=False)<=HTO_num].T  # filter unused HTOs

def index_mapping(indexs, index_map, index_origen, index_destination):
	# indexs as pd.series
	index_map.index = index_map[index_origen]
	indexs_mapped = indexs.apply(lambda x: index_map.loc[x.split('_')[0],index_destination]+'_'+x.split('_')[1]+'_'+x.split('_')[2]).values
	return indexs_mapped

# Filter based on cell barcodes
cell_barcode = pd.read_csv(data+prefix+'.edrops.matrix.hist.txt', sep='\t', header=None, index_col=0)  # Read in cell barcode
index_map = pd.read_csv(index_map_file, sep='\t', dtype='str')  # Read in RNA-seq and HTO sequencing index mapping table
fcts_matrix.index = index_mapping(fcts_matrix.index.to_series(), index_map, 'HTO','RNA-seq')  # convert HTO index to RNA-seq index
hashtag = fcts_matrix.reindex(cell_barcode.index, fill_value=0)  # filter unvalid cell barcodes

# Demultiplex HTOs using hashsolo
def hashsolo(ht_df):
	ht_adata = anndata.AnnData(obs = ht_df)  # create AnnData object

	#Run HashSolo via ScanPy
	htos = ht_df.columns
	sc.external.pp.hashsolo(ht_adata, htos)
	counts = ht_adata.obs['Classification'].value_counts()
	print(counts)  # print counts of HTOs
	
	# Generate output dataframe
	ht_adata.obs.to_csv(data+prefix+'_HTO.demultiplexed.txt', sep='\t')

hashsolo(hashtag)



