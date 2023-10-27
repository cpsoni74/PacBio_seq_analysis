import Bio.SeqIO
from Bio import SeqIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
import io
import random
import pysam
import time 
import math
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import copy
import time

def gen_mut_table(gbFile):
    bamfile = pysam.AlignmentFile('tempReadssorted.bam', 'rb', index_filename = 'tempReadssorted.bam.bai')
    startTime = time.time()
    mutationDict = {}
    MaPylRSDNA = 'ATGACGGTGAAATACACTGACGCACAGATCCAGCGCCTTCGCGAATATGGGAATGGCACGTATGAACAGAAAGTGTTCGAAGATTTGGCTTCGCGCGACGCAGCCTTTAGCAAAGAAATGAGTGTTGCCTCAACCGACAATGAGAAAAAAATTAAGGGCATGATTGCCAACCCGTCACGTCATGGACTTACGCAACTTATGAACGACATTGCCGACGCATTAGTCGCTGAGGGATTTATCGAGGTCCGCACGCCAATCTTTATCTCAAAAGACGCGCTTGCCCGTATGACGATTACAGAAGACAAGCCCCTGTTCAAGCAAGTATTCTGGATCGACGAGAAGCGTGCCTTACGCCCAATGTTGGCTCCAAATNNKTATTCCGTTATGCGTGATTTGCGTGACCACACCGACGGCCCAGTGAAGATTTTCGAGATGGGGAGCTGTTTTCGCAAGGAAAGTCACAGTGGCATGCATTTGGAGGAGTTCACGATGCTGNNKCTTNNKGATATGGGACCGCGTGGTGATGCGACAGAGGTTTTAAAAAATTACATTAGTGTTGTGATGAAAGCAGCGGGATTGCCCGATTATGATTTAGTCCAGGAAGAGAGTGACGTCTACAAAGAAACTATCGATGTTGAGATTAACGGGCAAGAAGTATGTAGCGCTGCTGTCGGACCCCATTATCTGGATGCTGCCCATGATGTGCATGAACCTTGGTCTGGTGCTGGTTTCGGTTTGGAGCGCTTATTAACCATTCGTGAGAAATATTCCACAGTAAAGAAAGGGGGGGCAAGTATCTCGTACCTGAACGGTGCAAAAATTAATTAA'
    MaPylRSAA = 'MTVKYTDAQIQRLREYGNGTYEQKVFEDLASRDAAFSKEMSVASTDNEKKIKGMIANPSRHGLTQLMNDIADALVAEGFIEVRTPIFISKDALARMTITEDKPLFKQVFWIDEKRALRPMLAPN?YSVMRDLRDHTDGPVKIFEMGSCFRKESHSGMHLEEFTML?L?DMGPRGDATEVLKNYISVVMKAAGLPDYDLVQEESDVYKETIDVEINGQEVCSAAVGPHYLDAAHDVHEPWSGAGFGLERLLTIREKYSTVKKGGASISYLNGAKIN*'
    for i, extracted_reads in enumerate(bamfile):
        # Run on subset as demonstration
        #if i > 469:
        #    break
        print (i)
        #if i % 20000 == 307:
        #    print(str(100*i/len(barcodesToConsensus))[:5] + \
        #          '% finished in ' + str(time.time() - startTime)[:5] + ' seconds')
        aligned_pairs = extracted_reads.get_aligned_pairs(with_seq = True)
        mutations = []
        refPosNotNone = 0
        for i, (read_pos, ref_pos, ref_base) in enumerate(aligned_pairs):
            # Prevent reference position from being interpreted as None
            if ref_pos != None:
                refPosNotNone = ref_pos
            if ref_base is not None:
                ref_base = ref_base.upper()
            if read_pos is None:
                read_base = None
            else:
                read_base = extracted_reads.query_sequence[read_pos].upper()
            if (read_base != ref_base) and (ref_base != 'N'):
                mutations = mutations + [[ref_base, refPosNotNone, read_base]]
                
        # This makes a dataframe where each barcode is paired with a list of lists of mutations
        mutationDict[extracted_reads.qname] = [mutations]
    testMutDF = pd.DataFrame(mutationDict).transpose().reset_index().rename(columns = {'index' : 'Read_Name', 
                                               0 : 'Mutation_List'})
    testMutDF.to_pickle('mutationCollection060523again.pkl')

    # Retreive mutation information

    mutationDF = pd.read_pickle('mutationCollection060523again.pkl')
    mutationDF.to_csv('mutationDF.csv')


    # Write out mutation logic for each barcode - save as .csv

    gbRec = SeqIO.read(gbFile, "genbank")
    featuresDict = {'Read_name' : []}
    featureLocs = {}
    for feature in gbRec.features:
        featuresDict[feature.type] = []
        featureLocs[feature.type] = [int(feature.location.start), int(feature.location.end)]
        
    def findMutLocs(mutList):
        mutLocs = []
        for mutation in mutList:
            #print ("mutation: ", mutation)
            mutLocs.append(mutation[1])
        return mutLocs

    #Backbone mutations    
    def findBBmut(row):
        mutLocs = findMutLocs(row['Mutation_List'])
        bbMutFound = False
        for position in mutLocs:
            try:
                #This is heavily based on the order of the sequence
                if (position < featureLocs['barcode'][0]): 
                    bbMutFound =True                 
                if (position > featureLocs['barcode'][1]) and (position < featureLocs['MaPylRS'][0]):
                    bbMutFound = True
                if (position > featureLocs['MaPylRS'][1]):
                    bbMutFound = True
                
            except:
                pass
        return bbMutFound

    def findInsertions(row):
        insFound = False
        for mutation in row['Mutation_List']:
            # If there is no assigned reference base that indicates an insertion
            if mutation[0] == None:
                insFound = True
        return insFound


    def findDeletions(row):
        delFound = False
        for mutation in row['Mutation_List']:
            # If there is no assigned read base that indicates a deletion
            if mutation[2] == None:
                delFound = True
        return delFound



    def translate(seq): 
        table = { 
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
        } 
        protein = "" 
        if len(seq)%3 == 0: 
            for i in range(0, len(seq), 3): 
                codon = seq[i:i + 3]
                if codon in table:
                    try:   
                        protein += table[codon] 
                    except:
                        protein = "_"
        return protein 


    def aaRSMut(row):
        mutLocs = findMutLocs(row['Mutation_List'])
        mutLocsInFrame = []
        affectedCodons = []
        aaRSMutList = []
        originalCodon = 'XXX'
        mutantCodon = 'YYY'
        originalAA = 'X'
        mutPosition = -1
        mutantAA = 'Y'
        aaMut = 'X0Y'
        for i, position in enumerate(mutLocs):
            if ((position >= featureLocs['MaPylRS'][0]) and (position <= featureLocs['L125'][0])
                or (position >= featureLocs['L125'][1]) and (position <= featureLocs['N166'][0])
                or (position >= featureLocs['N166'][1]) and (position <= featureLocs['V168'][0])
                or (position >= featureLocs['V168'][1]) and (position <= featureLocs['MaPylRS'][1])): 
                aaRSMutList.append(row['Mutation_List'][i])
                # position-1 gets the position with a 0 index
                # dividing by 3 and taking the floor gets the codon position
                # The codon position already has a 0 index
                mutPosition = math.floor((position-featureLocs['MaPylRS'][0])/3)
                if mutPosition not in (affectedCodons):
                    affectedCodons.append(mutPosition)
        if len(affectedCodons) == 0:
            aaMut = 'WT'
            return aaMut, originalAA, (mutPosition + 1), mutantAA
        
        try:
            if len(affectedCodons) == 1:
                originalCodon = MaPylRSDNA[affectedCodons[0] * 3 : affectedCodons[0] * 3 + 3]
                mutCodonSeqList = list(originalCodon)
                for mutationNumber, mutation in enumerate(aaRSMutList):
                    withinCodonPosition = (mutation[1] - featureLocs['MaPylRS'][0])%3
                    mutCodonSeqList[withinCodonPosition] = mutation[2]
                    mutantCodon = ''.join(mutCodonSeqList)
            else:
                aaMut = 'Multiple_mutations'
                return aaMut, originalAA, (mutPosition + 1), mutantAA
            
        except:
            aaMut = 'Ambiguous_mutation'
            return aaMut, originalAA, (mutPosition + 1), mutantAA
        
        #if mutantCodon not in PROGRAMMED_MUTATION_CODONS:
        #    aaMut = 'Illegal_mutation'
        #    return aaMut, originalAA, (mutPosition + 1), mutantAA
        
        if (originalCodon != 'XXX') and (mutantCodon != 'YYY'):
            originalAA = translate(originalCodon)
            mutantAA = translate(mutantCodon)
            if originalAA == mutantAA:
                aaMut = mutantAA = 'Silent mutation'
            elif mutantAA == '_':
                aaMut = mutantAA = 'Nonsense mutation'
            else:
                # Add one to get out of zero-indexing
                aaMut = (originalAA + str(mutPosition + 1) + mutantAA)
        return aaMut, originalAA, (mutPosition + 1), mutantAA

       
    mutationDF['BackboneMut'] = mutationDF.apply(lambda row: findBBmut(row), axis = 1)
    mutationDF['InsertionsFound'] = mutationDF.apply(lambda row: findInsertions(row), axis = 1)
    mutationDF['DeletionsFound'] = mutationDF.apply(lambda row: findDeletions(row), axis = 1)
    mutationDF['MaPylRSCodonMut'], mutationDF['originalAA'], mutationDF['AApos'], mutationDF['mutAA'] =  zip(*mutationDF.apply(aaRSMut, axis = 1))

    mutationDF.to_csv("mutationpackage.csv")

    aaRSdf = pd.read_csv('aaRS_df.csv')

    # Filter away barcodes that break the rules

    #betterMutsDF = mutationDF[
    #                          
    #                          ~mutationDF['MaPylRSmut']
    #                         ]

    betterMutsDF = mutationDF[
        (mutationDF['MaPylRSCodonMut'] != 'Multiple_mutations') &
        (mutationDF['MaPylRSCodonMut'] != 'Ambiguous_mutation') &
        (mutationDF['MaPylRSCodonMut'] != 'Nonsense mutation')

    ]

    merged_df = pd.merge(aaRSdf, betterMutsDF, on='Read_Name')
    #print (merged_df)
    merged_df[['barcode', 'L125/N166/V168', 'MaPylRSCodonMut']].to_csv('lookupTable.csv')
    merged_df


    # What percent of reads are of the "good" barcodes?
    uniqueBCDF = pd.read_csv('barcodecount_qscore.csv')
    usableReads = uniqueBCDF[uniqueBCDF['barcode'].isin(merged_df['barcode'])]['count'].sum()
    totalReads = uniqueBCDF['count'].sum()
    print ('usableReads:', usableReads)
    print(str(100*usableReads/totalReads)[:5] + '% of reads should be usable at t0')

    merged_df = pd.read_csv('lookupTable.csv')
    final_count = pd.DataFrame(merged_df['barcode'].value_counts())
    final_count.reset_index(inplace = True)
    final_count.rename(columns={'': 'index', 'index': 'barcode', 'barcode': 'count'}, inplace=True)

    active_site = merged_df.groupby('barcode')['L125/N166/V168'].apply(list).reset_index()
    Pyl_mut = merged_df.groupby('barcode')['L125/N166/V168'].apply(list).reset_index()

    barcoded_aaRS = pd.merge(final_count, active_site, on='barcode')
    barcoded_aaRS = pd.merge(final_count, Pyl_mut, on='barcode')
    barcoded_aaRS.to_csv("grouped by barcoded_aaRS.csv")

