import Bio.SeqIO
from Bio import SeqIO
from Bio import SeqFeature
import copy
""" 
---------------------------------------------------------------------------------
Load in .gb file and use it to find relevant features and locations with location\
---------------------------------------------------------------------------------
"""
def change_feature_type(template, i, strtloc, endloc, new_type):

    gbFile =  template
    gbRec = SeqIO.read(gbFile, "genbank")
    #print (featuresDict)
    # backbone is a join here, needs to be considered as multiple regions
    # Positions here are indexed to 1 so we will need to adjust for that too
    index = i

    feature_to_change = copy.deepcopy(gbRec.features[index])
    new_feature_type = SeqFeature.FeatureLocation(strtloc, endloc, 1) #create a new feature location object

    feature_to_change.type = new_type #change old feature location

    del gbRec.features[index] #remove changed feature

    gbRec.features.append(feature_to_change) #add identical feature with new location

    gbRec.features = sorted(gbRec.features, key = lambda feature: feature.location.start) # if you want them sorted by the start of the location, which is the usual case
    print (feature_to_change.type, "changed")
    SeqIO.write(gbRec, 'edited_template.gb', 'gb')
    edited_template = 'edited_template.gb'
    return edited_template
"""

#changing misc_feature to barcode, 2nd feature
index = 1

feature_to_change = copy.deepcopy(gbRec.features[index])
new_feature_type = SeqFeature.FeatureLocation(131, 140, 1) #create a new feature location object

feature_to_change.type = "barcode" #change old feature location

del gbRec.features[index] #remove changed feature

gbRec.features.append(feature_to_change) #add identical feature with new location

gbRec.features = sorted(gbRec.features, key = lambda feature: feature.location.start) # if you want them sorted by the start of the location, which is the usual case

print (feature_to_change)"""






