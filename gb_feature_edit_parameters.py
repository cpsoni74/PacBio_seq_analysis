#calling edit gb file with parameters to change feature type

import Bio.SeqIO
from Bio import SeqIO
from Bio import SeqFeature
import copy
from edit_gb_file_Feature import change_feature_type

#def change_feature_type(template, i, strtloc, endloc, new_type):
#template is gb file, i is the index of the feature, strtloc is the start of the feature as mentioned in gb file, endloc is end, new_type is the new annotation/new feature type

def annotating_gb (template):

	#10N barcode NNK library
	edit_1 = change_feature_type(template, 1, 131, 140, 'barcode')

	#L125 
	edit_2 = change_feature_type(edit_1, 6, 690, 692, 'L125')

	#N166 
	edit_3 = change_feature_type(edit_2, 7, 813, 815, 'N166')

	#V168 
	edit_4 = change_feature_type(edit_3, 8, 819, 821, 'V168')

	#MaPylRS
	edit_5 = change_feature_type(edit_4, 4, 318, 1145, 'MaPylRS')

	return edit_1

