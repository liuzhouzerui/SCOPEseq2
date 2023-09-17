#! /usr/bin/python
import gzip
import sys
import io

def clipper(read2_fastq,read2_clipped_fastq):
	i=0
	polya = 'AAAAAAAA'
	ct = 0
	cta = 0
	new=0
	with gzip.open(read2_clipped_fastq, 'wb') as g:
		with io.BufferedReader(gzip.open(read2_fastq,'rb')) as f:
			for line in f:
				dline = line.decode()
				if i == 1:
					pa = dline.find(polya)
					if pa!=-1:
						cta+=1
						if pa > 0:
							newline = dline[0:pa]+'\n'
						else:
							newline = 'NNN\n'
							pa=3
						newline = newline.encode()
						new = 1
						g.write(newline)
					else:
						g.write(line)
				elif i == 3:
					if new == 1:
						newline = dline[0:pa]+'\n'
						newline = newline.encode()
						new = 0
						g.write(newline)
					else:
						g.write(line)	
				else:
					g.write(line)
				i+=1
				if i == 4:
					i = 0
					new = 0
					ct+=1
	return cta,ct

