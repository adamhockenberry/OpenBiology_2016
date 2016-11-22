'''
File: adamAnalysis.py
Author: Adam Hockenberry
Description: 
'''

##Imports

#For biological sequence analysis
from Bio import SeqIO
from Bio.Seq import Seq

#Numerical methods
import numpy as np
from scipy import stats

#Plotting library
#import matplotlib as mpl
from matplotlib import pyplot as plt


def load_data(forward_riboprof_wig, reverse_riboprof_wig, forward_rnaseq_wig, reverse_rnaseq_wig, genome_genbank, simple_wig=True):
    profilingDict, fivePrimeProfile, sequenceDict, fivePrimeSeqDict, threePrimeSeqDict =\
                get_profiling_dict(forward_riboprof_wig, reverse_riboprof_wig, genome_genbank, simple_wig)
    MRNA, fivePrimeMRNA, sequenceDict, fivePrimeSeqDict, threePrimeSeqDict =\
                get_profiling_dict(forward_rnaseq_wig, reverse_rnaseq_wig, genome_genbank, simple_wig)
    return profilingDict, fivePrimeProfile, sequenceDict, fivePrimeSeqDict, MRNA, fivePrimeMRNA

def get_profiling_dict(forwardFileName, reverseFileName, genomeLocation, simple_wig=True):
    occupancyDictFwd = {}
    f = open(forwardFileName, 'r').readlines()
    if simple_wig:
        for line in f:
            tempLine = line.strip('\r\n').split('\t')
            occupancyDictFwd[int(tempLine[0])] = float(tempLine[1])
    else:
        counter = 1
        for line in f[2:]:
            tempLine = line.strip('\r\n')
            occupancyDictFwd[counter] = float(tempLine)
            counter += 1
    
    occupancyDictRev = {}
    f = open(reverseFileName, 'r').readlines()
    if simple_wig: 
        for line in f:
            tempLine = line.strip('\r\n').split('\t')
            occupancyDictRev[int(tempLine[0])] = float(tempLine[1])
    else:
        counter = 1
        for line in f[2:]:
            tempLine = line.strip('\r\n')
            occupancyDictRev[counter] = float(tempLine)
            counter += 1

    genome = list(SeqIO.parse(genomeLocation, 'gb'))[0]
    seqDict = {}
    fivePrimeSeqDict = {}
    threePrimeSeqDict = {}
    readsDict = {}
    fivePrimeReadsDict = {}
    for feature in genome.features:
        if feature.type == 'CDS':
            start = feature.location.start
            stop = feature.location.end
            locus = feature.qualifiers['locus_tag'][0]
            if feature.strand == 1:
                seqDict[locus] = str(genome.seq)[start:stop]
                fivePrimeSeqDict[locus] = str(genome.seq)[start-50:start]
                threePrimeSeqDict[locus] = str(genome.seq)[stop:stop+50]
                readsDict[locus] = get_occupancy(start, stop, occupancyDictFwd)       
                fivePrimeReadsDict[locus] = get_occupancy(start-50, start, occupancyDictFwd)       
            elif feature.strand == -1:
                seqDict[locus] = str(Seq(str(genome.seq)[start:stop]).reverse_complement())
                fivePrimeSeqDict[locus] = str(Seq(str(genome.seq)[stop:stop+50]).reverse_complement())
                threePrimeSeqDict[locus] = str(Seq(str(genome.seq)[start-50:start]).reverse_complement())
                readsDict[locus] = get_occupancy(start, stop, occupancyDictRev, reverse=True)       
                fivePrimeReadsDict[locus] = get_occupancy(stop, stop+50, occupancyDictRev, reverse=True)       
    return readsDict, fivePrimeReadsDict, seqDict, fivePrimeSeqDict, threePrimeSeqDict      

def get_occupancy(start, stop, occupancyDict, reverse=False):
    listy = []
    for position in range(start, stop):
        try:
            listy.append(occupancyDict[position])
        except KeyError:
            listy.append(0)
    if reverse:
        listy = listy[::-1]
    return listy


def example_profile(geneName, profilingDict, mrnaDict):
    fig = plt.figure(figsize=(9,6))
    ax1 = fig.add_subplot(211, frame_on=True, axisbg='white')
    ax1.plot(profilingDict[geneName])
    ax1.set_xlim(0, len(profilingDict[geneName]))
    ax1.set_xticklabels('')
    
    ax2 = fig.add_subplot(212, frame_on=True, axisbg='white')
    ax2.plot(mrnaDict[geneName])
    ax2.set_xlim(0, len(mrnaDict[geneName]))

    ax1.tick_params(labelsize=24)
    ax2.tick_params(labelsize=24)
    #plt.tight_layout()
    return


def get_efficiencies(profilingDict, mrnaDict, length_requirement, coverage_requirement, end_bases_to_ignore):
    total_prof = np.sum([inner for outer in profilingDict.values() for inner in outer])
    total_mrna = np.sum([inner for outer in mrnaDict.values() for inner in outer])
    transEffDict = {}
    mrnaExpDict = {}
    for gene in profilingDict.keys():
        if len(profilingDict[gene]) >= length_requirement:
            mrna_gene = mrnaDict[gene]
            if end_bases_to_ignore > 0:
                profiling_gene = profilingDict[gene][end_bases_to_ignore:-end_bases_to_ignore]
            else:
                profiling_gene = profilingDict[gene]
            
            if max(profilingDict[gene][:10]) > 0\
                    and np.percentile(mrna_gene, coverage_requirement) > 0\
                    and np.percentile(profiling_gene, coverage_requirement) > 0:
                        cds_rpkm = (sum(profiling_gene)*1000000000) / (len(profiling_gene)*total_prof)
                        mrna_rpkm = (sum(mrna_gene)*1000000000) / (len(mrna_gene)*total_mrna)
                        transEffDict[gene] = cds_rpkm / mrna_rpkm
                        mrnaExpDict[gene] = mrna_rpkm
    return transEffDict, mrnaExpDict

def normalize_dict_values(dicty):
    maximum = float(max(list(dicty.values())))
    for gene in dicty:
        dicty[gene] = dicty[gene] / maximum
    return dicty

def trans_eff_hist(transEffDictA, transEffDictB, transEffDictC, labels):
    #Regular scale
    bins=25
    y_A, binEdges_A = np.histogram(list(transEffDictA.values()), bins=bins, normed=True)
    bincenters_A = 0.5*(binEdges_A[1:]+binEdges_A[:-1])
    
    y_B, binEdges_B = np.histogram(list(transEffDictB.values()), bins=bins, normed=True)
    bincenters_B = 0.5*(binEdges_B[1:]+binEdges_B[:-1])
    
    y_C, binEdges_C = np.histogram(list(transEffDictC.values()), bins=bins, normed=True)
    bincenters_C = 0.5*(binEdges_C[1:]+binEdges_C[:-1])
    
    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(121)
    ax1.plot(bincenters_A, y_A, label='{} (n={})'.format(labels[0], len(transEffDictA.keys())), linewidth=3)
    ax1.plot(bincenters_B, y_B, label='{} (n={})'.format(labels[1], len(transEffDictB.keys())), linewidth=3)
    ax1.plot(bincenters_C, y_C, label='{} (n={})'.format(labels[2], len(transEffDictC.keys())), linewidth=3)
    ax1.autoscale()
    legend = ax1.legend(fontsize=16, loc='best').get_frame().set_linewidth(0.0)
    
    #Log scale
    bins=200
    y_A, binEdges_A = np.histogram(list(transEffDictA.values()), bins=bins, normed=True)
    bincenters_A = 0.5*(binEdges_A[1:]+binEdges_A[:-1])
    
    y_B, binEdges_B = np.histogram(list(transEffDictB.values()), bins=bins, normed=True)
    bincenters_B = 0.5*(binEdges_B[1:]+binEdges_B[:-1])
    
    y_C, binEdges_C = np.histogram(list(transEffDictC.values()), bins=bins, normed=True)
    bincenters_C = 0.5*(binEdges_C[1:]+binEdges_C[:-1])
    
    ax2 = fig.add_subplot(122)
    ax2.plot(bincenters_A, y_A, label=labels[0], linewidth=3)
    ax2.plot(bincenters_B, y_B, label=labels[1], linewidth=3)
    ax2.plot(bincenters_C, y_C, label=labels[2], linewidth=3)
    ax2.autoscale()
    ax2.set_xscale('log')
    plt.tight_layout()
    
    return
