import Bio.SeqIO
from Bio import SeqIO
from Bio import SeqFeature
import copy
""" 
---------------------------------------------------------------------------------
Load in .gb file and use it to find relevant features and locations with location\
---------------------------------------------------------------------------------
"""
def read_template(template):

    gbFile =  template
    gbRec = SeqIO.read(gbFile, "genbank")
    #print ("gbREc")
    #print (gbRec) 
    featuresDict = {'Read_name' : []}
    featureLocs = {}
    #print (featuresDict)
    # backbone is a join here, needs to be considered as multiple regions
    # Positions here are indexed to 1 so we will need to adjust for that too
    # Also pay attention to which features are reversed
    # I made sure to include stop codons in the two genes
   
    for feature in gbRec.features:
        featuresDict[feature.type] = []
        featureLocs[feature.type] = [int(feature.location.start), int(feature.location.end)]
        
    #print (featureLocs)
    #print (featuresDict)
    #print (gbRec.features)
    #print (len(gbRec))
    #print (featureLocs)
    return (gbRec, featuresDict, featureLocs)
 