#For biological sequence analysis
from Bio import SeqIO
from Bio.Seq import Seq

#Numerical methods
import numpy as np
from scipy import stats
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import summary_table

#Useful libraries
import random
import json

#Plotting library and associated settings to make pretty plots
from matplotlib import pyplot as plt

#For writing results
import datetime
import os

def get_seq_dicts(genomeLocation):
    genome = list(SeqIO.parse(genomeLocation, 'gb'))[0]
    seqDict = {}
    fivePrimeSeqDict = {}
    threePrimeSeqDict = {}
    for feature in genome.features:
        if feature.type == 'CDS':
            start = feature.location.start
            stop = feature.location.end
            locus = feature.qualifiers['locus_tag'][0]
            if feature.strand == 1:
                seqDict[locus] = str(genome.seq)[start:stop]
                fivePrimeSeqDict[locus] = str(genome.seq)[start-50:start]
                threePrimeSeqDict[locus] = str(genome.seq)[stop:stop+50]
            elif feature.strand == -1:
                seqDict[locus] = str(Seq(str(genome.seq)[start:stop]).reverse_complement())
                fivePrimeSeqDict[locus] = str(Seq(str(genome.seq)[stop:stop+50]).reverse_complement())
                threePrimeSeqDict[locus] = str(Seq(str(genome.seq)[start-50:start]).reverse_complement())
    return seqDict, fivePrimeSeqDict, threePrimeSeqDict

def get_binding_dict(five_prime_seq_dict, seq_dict, energyReference, asd_seq):
    temp_binding_dict = {}
    core_loc = asd_seq.index('CCUCC')
    for gene in five_prime_seq_dict.keys():
        temp_binding_dict[gene] = {}
        tempSeq = (five_prime_seq_dict[gene][20:]).replace('T', 'U')
        for event in range(len(tempSeq)-len(asd_seq)):
            temp_binding_value = energyReference[tempSeq[event:event+len(asd_seq)]]#NOTE: check this out later
            temp_binding_dict[gene][event-30+len(asd_seq)] = temp_binding_value
    return temp_binding_dict

def get_pandas_data_frame(trans_eff_dict, temp_binding_dict, spacing):
    xvals = []
    yvals = []
    names = [gene for gene in trans_eff_dict]
    a = [(trans_eff_dict[gene]) for gene in trans_eff_dict]
    b = [temp_binding_dict[gene][spacing] for gene in trans_eff_dict]

    table = list(zip(names,a,b))
    table.sort(key=lambda x: x[2])
    sorted_sdness = [z for x,y,z in table]
    headers = ['gene_names', 'teff', 'sdness']
    df = pd.DataFrame(table, columns=headers)
    df = df.set_index('gene_names')
    x_plotting = pd.DataFrame({'sdness': np.linspace(df.sdness.min(), df.sdness.max(), 100)})
    return df, sorted_sdness, x_plotting

def plot_individual(trans_eff_dict, five_prime_seq_dict, seq_dict, energy_file, spacing):
    asd_seq = energy_file.split('energyRef-')[-1].strip('.txt')
    with open(energy_file, 'r') as infile:
        energyReference = json.load(infile)
        temp_binding_dict = get_binding_dict(five_prime_seq_dict, seq_dict, energyReference, asd_seq)
        df, sorted_sdness, x_plotting = get_pandas_data_frame(trans_eff_dict, temp_binding_dict, spacing)
        
        ###Polynomial models
        poly_1 = smf.ols(formula='teff ~ sdness', data=df).fit()
        poly_2 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0)', data=df).fit()
        poly_3 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0) + I(sdness ** 3.0)', data=df).fit()
        poly_4 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0) + I(sdness ** 3.0) + I(sdness **4.0)', data=df).fit()

        poly_5 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0) + I(sdness ** 3.0) + I(sdness **4.0) + I(sdness ** 5.0)', data=df).fit()
        
        ###For Confidence Intervals
        st, data, ss2 = summary_table(poly_3, alpha=0.05)
        fitted_values_SD = data[:,2]
        predict_mean_se = data[:,3]
        predict_mean_ci_low_3, predict_mean_ci_up_3 = data[:, 4:6].T
        st, data, ss2 = summary_table(poly_1, alpha=0.05)
        fitted_values_SD = data[:,2]
        predict_mean_se = data[:,3]
        predict_mean_ci_low_1, predict_mean_ci_up_1 = data[:, 4:6].T

        pval_f_1st = poly_1.f_pvalue
        pval_f_1st = '%.2e' % pval_f_1st
        
        pval_f_3rd = poly_3.f_pvalue
        pval_f_3rd = '%.2e' % pval_f_3rd

        fig = plt.figure(figsize=(8,6))
        ax1 = fig.add_subplot(111)
        ax1.set_color_cycle(plt.cm.Dark2(np.linspace(0, 1, 2)))
        ax1.plot(df.sdness, df.teff, 'ko', alpha=0.2, label='')
        plt.xlabel('Binding strength (kcal/mol)', fontsize=24)
        #plt.ylabel('Residuals (log(RTE))', fontsize=24)
        plt.ylabel('log(RTE)', fontsize=24)
        plt.plot(x_plotting.sdness, poly_1.predict(x_plotting),\
                'b-', label='1st order ($R^{{2}}_{{adj}}$={:.3f}, $p$={})'.format(poly_1.rsquared_adj, pval_f_1st), alpha=0.9, linewidth=4)
        ax1.fill_between(sorted_sdness, predict_mean_ci_up_1, predict_mean_ci_low_1, facecolor='b', alpha=0.5)
        plt.plot(x_plotting.sdness, poly_3.predict(x_plotting),\
                'r-', label='3rd order ($R^{{2}}_{{adj}}$={:.3f}, $p$={})'.format(poly_3.rsquared_adj, pval_f_3rd), alpha=0.9, linewidth=4)
        ax1.fill_between(sorted_sdness, predict_mean_ci_up_3, predict_mean_ci_low_3, facecolor='r', alpha=0.5)
        plt.tick_params(labelsize=24)
        plt.legend(loc=4, fontsize=24).get_frame().set_linewidth(0.0)
        plt.tight_layout()
    too_strong = (poly_3.predict(x_plotting)[0])
    too_weak = (poly_3.predict(x_plotting)[-1])
    just_right = (max(poly_3.predict(x_plotting)))
    print (poly_1.summary())
    print (poly_3.summary())
    
    bin_edges = stats.mstats.mquantiles(df.sdness, [0, 0.2, 0.4, 0.6, 0.8, 1])

    bar_dicty = {}
    for j in range(len(bin_edges)-1):
        bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])] = []
    for i in df.index:
        for j in range(len(bin_edges)-1):
            if j != len(bin_edges)-2:
                if df.loc[i]['sdness'] >= bin_edges[j] and df.loc[i]['sdness'] < bin_edges[j+1]:
                    bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])].append(df.loc[i]['teff'])
                    break
            else:
                if df.loc[i]['sdness'] >= bin_edges[j]:
                    bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])].append(df.loc[i]['teff'])
                    break
    
    plot_list = []
    labels = ['']
    ns = []
    for j in range(len(bin_edges)-1):
        plot_list.append(bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])])
        n = len(bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])])
        labels.append('{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1]))
        ns.append(ns)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.axhline(0, c='k')
    ax.bar(range(len(bar_dicty)), [np.mean(i) for i in plot_list], yerr=[np.std(i)/len(i) for i in plot_list],\
            color=['0.2', '0.35', '0.50', '0.65', '0.8'], ecolor='k', align='center', zorder=3)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    plt.xlabel('Binding strength quintiles (kcal/mol)', fontsize=24)
    plt.tight_layout()
    
    xvals = [1, 2, 3, 4, 5]
    y_aic = [poly_1.aic, poly_2.aic, poly_3.aic, poly_4.aic, poly_5.aic]
    y_rsquared_adj = [poly_1.rsquared_adj, poly_2.rsquared_adj, poly_3.rsquared_adj, poly_4.rsquared_adj, poly_5.rsquared_adj]
    
    fig = plt.figure(figsize=(8,6))
    ax1 = fig.add_subplot(211)
    ax1.plot(xvals, y_rsquared_adj, marker='s', markersize=10)
    ax1.set_ylim(min(y_rsquared_adj)-(max(y_rsquared_adj)-min(y_rsquared_adj))*0.1,\
            max(y_rsquared_adj)+(max(y_rsquared_adj)-min(y_rsquared_adj))*0.1)
    ax1.set_xticklabels('')

    ax2 = fig.add_subplot(212)
    ax2.plot(xvals, y_aic, marker='s', markersize=10)
    ax2.set_ylim(min(y_aic)-(max(y_aic)-min(y_aic))*0.1, max(y_aic)+(max(y_aic)-min(y_aic))*0.1)
    plt.tight_layout()
    return poly_1, poly_3, percents, df

def compare_first_and_third_degree(trans_eff_dict, five_prime_seq_dict, seq_dict, energy_file, spacing_list):
    asd_seq = energy_file.split('energyRef-')[-1].strip('.txt')
    poly_1_list = []
    poly_3_list = []
    spacing_plot = []
    for spacing in spacing_list:
        spacing = spacing*-1
        asd_seq = energy_file.split('energyRef-')[-1].strip('.txt')
        with open(energy_file, 'r') as infile:
            energyReference = json.load(infile)
            temp_binding_dict = get_binding_dict(five_prime_seq_dict, seq_dict, energyReference, asd_seq)
            df, sorted_sdness, x_plotting = get_pandas_data_frame(trans_eff_dict, temp_binding_dict, spacing)

            poly_1 = smf.ols(formula='teff ~ sdness', data=df).fit()
            poly_3 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0) + I(sdness ** 3.0)', data=df).fit()
            poly_1_list.append(poly_1.rsquared_adj)
            poly_3_list.append(poly_3.rsquared_adj)
            spacing_plot.append(spacing)
    
    fig = plt.figure(figsize=(8,6))
    ax1 = fig.add_subplot(111)
    ax1.plot(spacing_plot, poly_1_list, 'b-', label='1st order', linewidth=4, alpha=0.6)
    ax1.plot(spacing_plot, poly_3_list, 'r-', label='3rd order', linewidth=4, alpha=0.6)
    ax1.plot(spacing_plot, poly_1_list, 'bs', markersize=10)
    ax1.plot(spacing_plot, poly_3_list, 'rs', markersize=10)
    ax1.axhline(0.0, linewidth=2, c='k')
    plt.xlabel('Distance to start codon', fontsize=24)
    plt.ylabel('$R_{adj}^{2}$', fontsize=24)
    plt.legend(['1st order', '3rd order'], loc='best', fontsize=24).get_frame().set_linewidth(0.0)
    plt.tick_params(labelsize=24)
    plt.tight_layout()
    return

def compare_different_asds(trans_eff_dict, five_prime_seq_dict, seq_dict, energy_file_list, spacing_list):
    asd_seqs = [] 
    poly_3_meta = []
    for energy_file in energy_file_list:
        print(energy_file)
        asd_seq = energy_file.split('energyRef-')[-1].strip('.txt')
        with open(energy_file, 'r') as infile:
            energyReference = json.load(infile)
            temp_binding_dict = get_binding_dict(five_prime_seq_dict, seq_dict, energyReference, asd_seq)
        poly_3_list = []
        spacing_plot = []
        for spacing in spacing_list:
            spacing = spacing*-1
            asd_seq = energy_file.split('energyRef-')[-1].strip('.txt')
            df, sorted_sdness, x_plotting = get_pandas_data_frame(trans_eff_dict, temp_binding_dict, spacing)
            poly_3 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0) + I(sdness ** 3.0)', data=df).fit()
            poly_3_list.append(poly_3.rsquared_adj)
            spacing_plot.append(spacing)
                
        poly_3_meta.append(poly_3_list[::-1])
        asd_seqs.append(asd_seq)
    fig = plt.figure(figsize=(8,6))
    ax1 = fig.add_subplot(111)
    cax = ax1.matshow(poly_3_meta, cmap='BuPu')
    ax1.set_yticklabels('')
    cbar = fig.colorbar(cax)
    plt.tight_layout()
    return

def plot_taniguchi_data(first_order_r2, third_order_r2, percent_increase_list, x_labels):
    spacing = np.arange(len(first_order_r2))
    width = 0.35
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.bar(spacing, first_order_r2, width, color='blue', alpha=0.6, label='First order', align='center', zorder=3)
    ax.bar(spacing+width, third_order_r2, width, color='red', alpha=0.6, label='Third order', align='center', zorder=3)
    ax.legend(loc='best', fontsize=24).get_frame().set_linewidth(0.0)
    ax.set_xticklabels(x_labels, rotation=45)
    ax.set_ylabel('R^{2}_{adj}', size=32)
    ax.set_xlabel('Signal/Noise quality thresholds', size=32)
    plt.tight_layout()
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.bar(spacing, np.array(percent_increase_list)*100, color=['0.8', '0.65', '0.5', '0.35', '0.2'], align='center', zorder=3)
    ax.set_xticklabels(x_labels, rotation=45)
    ax.set_ylabel('\% increase in mean,\n 1st to 4th quintiles', size=32)
    ax.set_xlabel('Signal/Noise quality thresholds', size=32)
    plt.tight_layout()
    return

def plot_individual_kosuri(exp_dicty, gene_binding_dict, asd_seq, spacing):
    teff_list = []
    sd_list = []
    for i in exp_dicty.keys():
        if i != 'DeadRBS':
            #if i not in ['GSGV_RBS', 'GSG_RBS']:
            #if i not in ['GSG_RBS']:
            #if i not in ['GSGV_RBS']:
                teff_list.append(exp_dicty[i])
                sd_list.append(gene_binding_dict[i])
    
    table = list(zip(teff_list,sd_list))
    table.sort(key=lambda x: x[1])
    sorted_sdness = [y for x,y in table]
    headers = ['teff', 'sdness']
    df = pd.DataFrame(table, columns=headers)
    x_plotting = pd.DataFrame({'sdness': np.linspace(df.sdness.min(), df.sdness.max(), 100)})

    ###Polynomial models
    poly_1 = smf.ols(formula='teff ~ sdness', data=df).fit()
    poly_2 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0)', data=df).fit()
    poly_3 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0) + I(sdness ** 3.0)', data=df).fit()
    poly_4 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0) + I(sdness ** 3.0) + I(sdness **4.0)', data=df).fit()
    poly_5 = smf.ols(formula='teff ~ sdness + I(sdness ** 2.0) + I(sdness ** 3.0) + I(sdness **4.0) + I(sdness ** 5.0)', data=df).fit()
    
    ###For Confidence Intervals
    st, data, ss2 = summary_table(poly_3, alpha=0.05)
    fitted_values_SD = data[:,2]
    predict_mean_se = data[:,3]
    predict_mean_ci_low_3, predict_mean_ci_up_3 = data[:, 4:6].T
    st, data, ss2 = summary_table(poly_1, alpha=0.05)
    fitted_values_SD = data[:,2]
    predict_mean_se = data[:,3]
    predict_mean_ci_low_1, predict_mean_ci_up_1 = data[:, 4:6].T

    pval_f_1st = poly_1.f_pvalue
    pval_f_1st = '%.2e' % pval_f_1st
    pval_f_3rd = poly_3.f_pvalue
    pval_f_3rd = '%.2e' % pval_f_3rd

    fig = plt.figure(figsize=(8,6))
    ax1 = fig.add_subplot(111)
    ax1.set_color_cycle(plt.cm.Dark2(np.linspace(0, 1, 2)))
    ax1.plot(df.sdness, df.teff, 'ko', label='', alpha=0.2)
    plt.xlabel('Binding strength (kcal/mol)', fontsize=24)
    plt.ylabel('log(RTE)', fontsize=24)
    plt.plot(x_plotting.sdness, poly_1.predict(x_plotting),\
            'b-', label='1st order ($R^{{2}}_{{adj}}$={:.3f}, $p$={})'.format(poly_1.rsquared_adj, pval_f_1st), alpha=0.9, linewidth=4)
    ax1.fill_between(sorted_sdness, predict_mean_ci_up_1, predict_mean_ci_low_1, facecolor='b', alpha=0.5)
    plt.plot(x_plotting.sdness, poly_3.predict(x_plotting),\
            'r-', label='3rd order ($R^{{2}}_{{adj}}$={:.3f}, $p$={})'.format(poly_3.rsquared_adj, pval_f_3rd), alpha=0.9, linewidth=4)
    ax1.fill_between(sorted_sdness, predict_mean_ci_up_3, predict_mean_ci_low_3, facecolor='r', alpha=0.5)
    
    plt.tick_params(labelsize=24)
    plt.legend(loc='best', fontsize=24).get_frame().set_linewidth(0.0)
    plt.tight_layout()
    
    too_strong = (poly_3.predict(x_plotting)[0])
    too_weak = (poly_3.predict(x_plotting)[-1])
    just_right = (max(poly_3.predict(x_plotting)))
    print (poly_1.summary())
    print (poly_3.summary())
    
    bin_edges = stats.mstats.mquantiles(df.sdness, [0, 0.2, 0.4, 0.6, 0.8, 1])
    #bin_edges = stats.mstats.mquantiles(df.sdness, [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

    bar_dicty = {}
    for j in range(len(bin_edges)-1):
        bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])] = []
    for i in range(len(df)):
        for j in range(len(bin_edges)-1):
            if j != len(bin_edges)-2:
                if df.loc[i]['sdness'] >= bin_edges[j] and df.loc[i]['sdness'] < bin_edges[j+1]:
                    bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])].append(df.loc[i]['teff'])
                    break
            else:
                if df.loc[i]['sdness'] >= bin_edges[j]:
                    bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])].append(df.loc[i]['teff'])
                    break
    
    plot_list = []
    labels = ['']
    ns = []
    for j in range(len(bin_edges)-1):
        plot_list.append(bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])])
        labels.append('{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1]))
        ns.append(len(bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[j], bin_edges[j+1])]))

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.axhline(0, c='k')
    ax.bar(range(len(bar_dicty)), [np.mean(i) for i in plot_list], yerr=[np.std(i)/len(i) for i in plot_list],\
            color=['0.2', '0.35', '0.50', '0.65', '0.8'], ecolor='k', align='center', zorder=3)
    ax.set_xticklabels(labels, rotation=45, ha='right')
    plt.xlabel('Binding strength quintiles (kcal/mol)', fontsize=24)
    plt.ylabel('log(RTE)', fontsize=24)
    plt.tight_layout()
    
    non_sd_bin = bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[4], bin_edges[5])]
    sd_bin = bar_dicty['{:.2f} to {:.2f}'.format(bin_edges[1], bin_edges[2])]
    percents = np.power(2, np.mean(sd_bin))/np.power(2, np.mean(non_sd_bin))-1
    print('#####################################################')
    print('Rank sum statistics between 2nd and 5th quintile:')
    print('Mean SD:', np.power(2, np.mean(sd_bin)))
    print('Mean non-SD:', np.power(2, np.mean(non_sd_bin)))
    print('Percent diff:', percents)
    print('Rank sums:', stats.ranksums(sd_bin, non_sd_bin))
    print('#####################################################')
    print('AICs probability that first order is better than third order:')
    print(np.power(2, (poly_3.aic-poly_1.aic)/2.))

    xvals = [1, 2, 3, 4, 5]
    y_aic = [poly_1.aic, poly_2.aic, poly_3.aic, poly_4.aic, poly_5.aic]
    y_rsquared_adj = [poly_1.rsquared_adj, poly_2.rsquared_adj, poly_3.rsquared_adj, poly_4.rsquared_adj, poly_5.rsquared_adj]
    
    fig = plt.figure(figsize=(8,6))
    ax1 = fig.add_subplot(211)
    ax1.plot(xvals, y_rsquared_adj, marker='s', markersize=10)
    ax1.set_ylim(min(y_rsquared_adj)-(max(y_rsquared_adj)-min(y_rsquared_adj))*0.1,\
            max(y_rsquared_adj)+(max(y_rsquared_adj)-min(y_rsquared_adj))*0.1)

    ax2 = fig.add_subplot(212)
    ax2.plot(xvals, y_aic, marker='s', markersize=10)
    ax2.set_ylim(min(y_aic)-(max(y_aic)-min(y_aic))*0.1, max(y_aic)+(max(y_aic)-min(y_aic))*0.1)
    plt.tight_layout()
    return



def teff_given_structure(transEffDict, structureDict):
    transEffList = [np.log2(transEffDict[i]) for i in transEffDict]
    structureList = [structureDict[i] for i in transEffDict]
    a, b, c, d, e = stats.linregress(structureList, transEffList)
    
    print(a,b,c,d,e, '***', c**2)


    correctedDict = {}
    for gene in transEffDict:
        correctedDict[gene] = np.log2(transEffDict[gene]) - (a*structureDict[gene] + b)

    correctedList = [correctedDict[i] for i in transEffDict]
    a, b, c, d, e = stats.linregress(structureList, correctedList)
    print(a,b,c,d,e)
    
    fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(121)
    ax1.plot(structureList, transEffList, color='steelblue', marker='o', linestyle='', alpha=0.5, markeredgewidth=0.2)
    labels = ax1.get_xticklabels()
    for label in labels:
        label.set_rotation(45)
    ax2 = fig.add_subplot(122)
    ax2.plot(structureList, correctedList, color='steelblue', marker='o', linestyle='', alpha=0.5, markeredgewidth=0.2)
    labels = ax2.get_xticklabels()
    for label in labels:
        label.set_rotation(45)
    plt.tight_layout()
    return correctedDict

def get_structure(sequenceDict, fivePrimeSeqDict, transEffDict):
    mfeStructureDict = {}
    for gene in transEffDict.keys():
        tempSeq = fivePrimeSeqDict[gene][20:] + sequenceDict[gene][:30]
        tempSeq = tempSeq.replace('T', 'U')
        write_temp_sequence(tempSeq, './')
        mfe, ensemblemfe, diversity = RNA_fold_wrapper('./')
        mfeStructureDict[gene] = mfe
    return mfeStructureDict

def write_temp_sequence(sequence, directory):
    temp = open(directory + 'temp.seq', 'w')
    temp.write('>' + 'temp' + '\n')
    temp.write(sequence)
    temp.close()


#################################################################################
#################################################################################
#Temporary rna folding stuff
#################################################################################
#################################################################################

def RNA_fold_wrapper(directory):
    f = os.popen("RNAfold -p < " + directory + "temp.seq")
    output = f.readlines()
    MFE = output[2]
    ensembleFE = output[3]
    diversity = output[-1][-15:]
    MFE = str(find_number(MFE))
    ensembleFE = str(find_number(ensembleFE))
    ensembleDiversity = str(find_number(diversity))
    return float(MFE), float(ensembleFE), float(ensembleDiversity)

def find_number(string):
    digits = 10
    temp = 'fail'
    for x in range(digits, 0, -1):
        for y in range(len(string)-x):
            try:
                temp = float(string[y:(y+x)])
            except ValueError:
                temp = 'fail'
            if type(temp).__name__ == 'float':
                break
        if type(temp).__name__ == 'float':
            break
    return temp
