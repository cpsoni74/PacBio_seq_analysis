#Pacbio Q-score analysis and barcode filtering

import pysam
import pandas as pd
import numpy as np 
from Bio import SeqIO

def qscore(temp_rec, bamfile, baifile):
	print (str(bamfile))
	BC = {'barcode' : []}
	QC = {'base' : [], 'phred' : []}
	QC2 = {'Read_Name': [], 'barcode' : [], 'phred' : []}

	featureLocs = {}
	for feature in temp_rec[0].features: 
		featureLocs[feature.type] = [int(feature.location.start) + 1, int(feature.location.end) + 1] 
	n = 1
	tot = 700000
	bclen = 10
	with pysam.AlignmentFile(bamfile, "rb", index_filename = baifile) as bam:
		BCReads = pysam.AlignmentFile("tempReads.bam", "wb", template = bam)
		for i, read in enumerate(bam.fetch()):
			barcode = ''
			QC['base'] = []
			QC['phred'] = []
			if i >tot:
				break
        	
			print ("_________________________________________- i:", i)

			if (read.flag == 0) | (read.flag == 16):
				#print ('hi')
				for alignedPair in read.get_aligned_pairs(with_seq = True):
					refPos = alignedPair[1]
					try:
						if ((refPos >= featureLocs['barcode'][0] - 1) & (refPos < featureLocs['barcode'][1] - 1)):
							if (read.query_qualities[refPos] >= 30):
								barcode = barcode + read.seq[alignedPair[0]]
								QC['base'].append(read.query_sequence[alignedPair[0]]) 
								QC['phred'].append(read.query_qualities[alignedPair[0]])							
					except:
						break
			else:
				pass

			if len(QC['base']) == bclen:
				n = n + 1
				QC2['Read_Name'].append(read.qname)
				QC2['barcode'].append(barcode)
				QC2['phred'].append(QC['phred'])
				BCReads.write(read)
			else:
				print ("low qual read for ", barcode, read.qname)
				pass
			print (barcode, len(barcode))
			BC['barcode'].append(barcode)
	print ('Number of incorrect  barcodes :',  i - n - 1) #i is 0
	print ('Correct barcodes : ', n)

	qual = pd.DataFrame(QC2)
	qual.to_csv('barcode_qscore.csv')

	qual = pd.read_csv('barcode_qscore.csv')
	uniqueBCDF = pd.DataFrame(qual['barcode'].value_counts())
	uniqueBCDF.reset_index(inplace = True)
	uniqueBCDF.rename(columns={'': 'index', 'index': 'barcode', 'barcode': 'count'}, inplace=True)

	BCDF_count = pd.DataFrame(uniqueBCDF)
	BCDF_count.to_csv('barcodecount_qscore.csv')
	BCDF_count

	BCReads.close()
