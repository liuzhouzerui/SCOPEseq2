# SCOPEseq2


## About
SCOPE-seq2 is an integrated single cell live imaging and RNA sequencing method. This pipeline is created for processing scRNA-seq/CITE-seq/cell hashing data generated with optically decodable barcoded beads in SCOPE-seq2.  
- [Paper at Scientific Reports (Liu *et al.* 2020)](https://doi.org/10.1038/s41598-020-76599-w)
- [Application to neurogenesis (Mizrak *et al.* 2020)](https://doi.org/10.1016/j.celrep.2020.107805)
- [Application to fluorescence guided surgery (Liu *et al.* 2022)](https://doi.org/10.1101/2022.12.17.520870)


## Installation dependnecies
- requires python3.6, bx-python, pysam
- STAR v2.7.0d
- pandas
- matplotlib
- optional: rpy2, DropletUtils (R) to run `emptyDrops`
- optional: anndata, scanpy to run `HashSolo`

## Usage
1. Processing illumina paired-end fastq files for scRNA-seq
```
python SCOPEseq2/SCOPEseq2_pipeline.py -bc1 SCOPEseq2/reference/S_barcode.csv -bc2 SCOPEseq2/reference/Q_barcode.csv -ind SCOPEseq2/reference/nextera_index.csv -ref STAR-REFERENCE-GENOME -gtf STAR-REFERENCE-GTF -t THREAD -r PREFIX -data OUTDIR -r1 READ1-FASTQ -r2 READ2-FASTQ
```
Output files

- `OURDIR/PREFIX.exon_address.filt.txt` A table for molecule counts

Optional: identify empty drops using [emptyDrops](https://doi.org/10.1186/s13059-019-1662-y)
```
python SCOPEseq2/edrops.py -f OURDIR/PREFIX.exon_address.filt.txt -r STAR-REFERENCE-GENOME/geneInfo.tab -p OURDIR/PREFIX -l 100
```
Output files

- `OURDIR/PREFIX.edrops.matrix.txt` Expression count matrix (gene x cell)


2. Processing illumina paired-end fastq files for cell hashing. This code can also be used for CITE-seq data to generate ADT counts.
```
python SCOPEseq2/SCOPEseq2_HTO_pipeline.py -bc1 SCOPEseq2/reference/S_barcode.csv -bc2 SCOPEseq2/reference/Q_barcode.csv -ind SCOPEseq2/reference/TruSeq_index.csv -b SCOPEseq2/reference/TotalseqA_HTO_human.csv -r PREFIX -data OUTDIR -r1 READ1-FASTQ -r2 READ2-FASTQ
```
Output files

- `OURDIR/PREFIX_HTO.citeseq_address.filt.txt` A table for HTO molecule counts

Optional: identify empty drops using [HashSolo](https://doi.org/10.1016/j.cels.2020.05.010)
```
python SCOPEseq2/HTO_demultiplex.py -p PREFIX -data OUTDIR -n NUMBER-OF-HTOs -ind INDEX-CONVERT.txt
```
where INDEX-CONVERT.txt is a whitespace delimited table with two columns `HTO RNA-seq`

Output files

- `OURDIR/PREFIX_HTO.demultiplexed.txt` A table for HTO sample / Doublet assignment
