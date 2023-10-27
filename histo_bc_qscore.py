#histogram for barcodes based on qscore

import matplotlib
import matplotlib.pyplot as plt
import copy
import time

# Histogram with integral that is # of reads
# Find how many barcodes have more than 1 read
def histo(final_BCDF, readsBCDF):

	readThreshold = 1
	usabilityMinimum = len(readsBCDF)
	print (usabilityMinimum)

	plt.gcf().set_size_inches((20, 5))  
	readBChistData = readsBCDF.merge(final_BCDF, left_on = 'barcode', right_on = 'barcode')
	binSet = range(0,50,1)
	plt.hist(readBChistData['count'], bins = binSet, histtype = 'step', color = 'black')
	plt.yscale('log')
	plt.ylabel('Reads of barcodes with that many reads')
	plt.xlabel('Reads per barcode')
	plt.axvline(x = readThreshold, linestyle = '--')
	plt.text(5, 2500, 'Usability threshold')
	plt.text(5, 500, 'Potentially usable barcodes: ' + str(usabilityMinimum))
	plt.savefig('barcodes_histogram_qscore.png', dpi = 600)
	#plt.show()

