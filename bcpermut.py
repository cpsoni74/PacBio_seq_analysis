#barcodes per mutation 
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def barcodes_per_mutation():

	merged_df = pd.read_csv('lookupTable.csv')
	mutation_counts = merged_df.groupby('L125/N166/V168')['barcode'].count().reset_index()
	mutation_counts.columns = ['mutation_count', 'barcode_count']
	#mutation_counts = mutation_counts.drop(index=mutation_counts.index[:71])
	mutation_counts = mutation_counts.reset_index(drop = True)
	mutation_counts.to_csv("aaRS vs barcode.csv")
	print (mutation_counts)

	# Sort the DataFrame based on the 'mutation_count' column
	mutation_counts = mutation_counts.sort_values('mutation_count')

	# Plot the bar plot
	plt.clf()
	binSet = range(0,50,1)
	plt.hist(mutation_counts['barcode_count'], bins = binSet, color = 'lightgray')
	plt.xlabel('Barcodes per mutation')
	plt.ylabel('Mutation')
	plt.title('Barcodes per Mutation Histogram')
	plt.savefig('aaRS_barcoded library.png', dpi = 600)
	#plt.show()
