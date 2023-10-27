#extracting aaRS 

import pysam 
import pandas as pd
import numpy as np 
from Bio import SeqIO
import re
from translate import translate

def extract_aaRS(temp_rec, qcutoff):
	final_lib = {'barcode': [], 'Read_Name': [], 'aaRS' : [], 'L125/N166/V168' : [],'phred' :[]}
	bamfile = 'tempReadssorted.bam'
	baifile = 'tempReadssorted.bam.bai'
	featureLocs = {}
	for feature in temp_rec[0].features: 
		featureLocs[feature.type] = [int(feature.location.start) + 1, int(feature.location.end) + 1] 
	
	bc_qsorted = pd.read_csv('barcode_qscore.csv')
	#bc_qsorted.update({'aaRS' : []})
	bam = pysam.AlignmentFile(bamfile, 'rb', index_filename = baifile)
	for i, read in enumerate(bam):
		print ('Round :', i)
		L125 = ''
		N166 = ''
		V168 = ''
		AA125 = ''
		AA166 = ''
		AA168 = ''
		aaRS_lib = []
		qc_lib = []
		barcode = ''
		AA = []
		x = 0
		
		#for barcode_read in bc_qsorted['Read_Name']:
		#	barcode_read_cleaned = re.sub(r'^\'|\'$', '', barcode_read.strip('[]'))  #removes [''] from the barcode_read because it is an element in the list
			#print (read.qname, barcode_read_cleaned)
		#	if barcode_read_cleaned == read.qname:
				#print (barcode_read_cleaned)
		dont = True
		for alignedPair in read.get_aligned_pairs(with_seq = True):
			refPos = alignedPair[1]
			#print ("refPos, alignedPair :", refPos, alignedPair)
			#print ("x =" , x)
			#if (i < 3034):
			try:
                # Subtract 1 for 0 index
                #print ("featureLocs")
                #print (featureLocs["barcode"])
                #print (refPos)
					#print ("refPos, alignedPair :", refPos, alignedPair)
				if ((refPos >= featureLocs['barcode'][0] - 1) and (refPos < featureLocs['barcode'][1] - 1)):
					barcode = barcode + read.seq[alignedPair[0]]
				if ((refPos >= featureLocs['L125'][0] - 1) and (refPos < featureLocs['L125'][1] - 1)):
					if (read.query_qualities[refPos] >= qcutoff):
						L125 = L125 + read.seq[alignedPair[0]]
						q125 = read.query_qualities[alignedPair[0]]	
						#print ("q125:" , q125)						
				if ((refPos >= featureLocs['N166'][0] - 1) and (refPos < featureLocs['N166'][1] - 1)):
					if (read.query_qualities[refPos] >= qcutoff):
						N166 = N166 + read.seq[alignedPair[0]]
						q166 = read.query_qualities[alignedPair[0]]	
						#print ("q166", q166)
				if ((refPos >= featureLocs['V168'][0] - 1) and (refPos < featureLocs['V168'][1] - 1)):
					if (read.query_qualities[refPos] >= qcutoff):
						V168 = V168 + read.seq[alignedPair[0]]
						q168 = read.query_qualities[alignedPair[0]]
						#print ("V168", V168)
					#print ("pass")	

			except:
				if (x < 3033):
					#basically the final base 3034 will be always be none, so it will enter the except loop regardless, don't worry about this. This is to prevent truncated sequences/low quality sequences to get in
					print ("fail", x)
					print(read.qname)
					dont = False
					break
			x = x + 1
		#print ("after inner loop")
		#print ("x =" , x)
			#print (L125)
		if dont:
			AA125 = translate(L125)
			AA166 = translate(N166)
			AA168 = translate(V168)
			if AA125 != '' and AA166 != '' and AA168 != '' and AA125 != '_' and AA166 != '_' and AA168 != '_':
				aaRS_lib.extend([L125, N166, V168])
				qc_lib.extend([q125, q166, q168])
				AA.extend([AA125, AA166, AA168])
				#print (dont, AA)
				final_lib['aaRS'].append(aaRS_lib)
				final_lib['barcode'].append(barcode)
				final_lib['Read_Name'].append(read.qname)
				final_lib['phred'].append(qc_lib)
				final_lib['L125/N166/V168'].append(AA)

	print (final_lib)
	aaRS_df = pd.DataFrame(final_lib)
	aaRS_df.to_csv('aaRS_df.csv')
	aaRS_df = pd.read_csv('aaRS_df.csv')
	count_df = pd.DataFrame(aaRS_df['barcode'].value_counts())
	count_df.reset_index(inplace = True)
	count_df.rename(columns={'': 'index', 'index': 'barcode', 'barcode': 'count'}, inplace=True)
	grouped_aaRSreads = aaRS_df.groupby('barcode')['Read_Name'].apply(list).reset_index()
	grouped_aaRS = aaRS_df.groupby('barcode')['aaRS'].apply(list).reset_index()
	qc_aaRS = aaRS_df.groupby('barcode')['phred'].apply(list).reset_index()
	AA_aaRS = aaRS_df.groupby('barcode')['L125/N166/V168'].apply(list).reset_index()

	aaRS_final = pd.merge(count_df, grouped_aaRSreads, on='barcode')
	aaRS_final = pd.merge(aaRS_final, grouped_aaRS, on='barcode')
	aaRS_final = pd.merge(aaRS_final, AA_aaRS, on='barcode')
	aaRS_final = pd.merge(aaRS_final, qc_aaRS, on='barcode')

	aaRS_final.to_csv('aaRS_count.csv', index=True)
 

