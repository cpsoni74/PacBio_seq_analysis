import pandas as pd
import matplotlib.pyplot as plt

def aaRS_per_barcode():
    merged_df = pd.read_csv('lookupTable.csv')

    # Group by 'barcode' and count unique values in 'L125/N166/V168'
    mutation_counts = merged_df.groupby('barcode')['L125/N166/V168'].nunique().reset_index()
    mutation_counts.columns = ['barcode', 'unique_mutation_count']


    # Count how many barcodes have a specific number of unique mutations
    counts = mutation_counts['unique_mutation_count'].value_counts().reset_index()
    counts.columns = ['unique_mutation_count', 'barcode_count']

    # Sort the counts by the number of unique mutations
    counts = counts.sort_values('unique_mutation_count')
    mutation_counts.to_csv("/Users/chintansoni/Desktop/MaPyl Barcoding library seq/PacBio input/Code/aaRSperbc.csv")

    # Plot the bar chart
    plt.clf()
    binSet = range(0, 10, 1)
    plt.hist(mutation_counts['unique_mutation_count'], bins=binSet, color='lightgray')
    plt.xlabel('Number of aaRS Mutations')
    plt.ylabel('Number of Barcodes')
    plt.title('Number of barcodes vs Number of aaRS mutations')
    plt.savefig('aaRS_perbc.png', dpi=600)
    #plt.show()

