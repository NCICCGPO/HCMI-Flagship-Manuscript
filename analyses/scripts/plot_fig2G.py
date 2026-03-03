import sys
import os
import random
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from pybedtools import BedTool
from itertools import product
from collections import Counter
from pybiomart import Server 
from scipy.stats import linregress
import h5py
import tables
from upsetplot import plot, UpSet
import upsetplot
import random





#Define a function to load the sample list
def load_sample_list():
    sample_list = pd.read_csv("seongmin_sample_map.csv")
    return sample_list






#Load the mutation signature cosine similarity
def load_mutsig_data():
    mutsig_data = pd.read_csv("mutation_sig_cosine_similarity_scores.csv", sep = ",")
    mutsig_data['unique_id'] = mutsig_data['left_sample'] + '_' + mutsig_data['right_sample']
    return mutsig_data




#Load in the consensus data from tim
def load_consensus_data():
    consensus_data = pd.read_csv("pairwise_loh.csv", sep = ",")
    consensus_data['unique_id'] = consensus_data['tumor_id'] + '_' + consensus_data['model_id']
    return consensus_data





#Load in the snv concordance data from tim
def load_snv_consensus_data():
    snv_consensus_data = pd.read_csv("SNV_concordance.txt", sep = "\t")
    snv_consensus_data.columns = ['unique_id', 'intersecting_snv', 'tumor_only_snv', 'model_only_snv']

    return snv_consensus_data





def construct_consensus_barplot(sample_list, mutsig_data, consensus_data, snv_consensus_data):
    joined_df = pd.merge(consensus_data, mutsig_data, on='unique_id', how='inner')
    joined_df = pd.merge(joined_df, snv_consensus_data, on='unique_id', how='inner')

    upset_plot_df = joined_df[['unique_id', 'loh_shared', 'model_consensus_purity', 'cosine_similarity', 'wgd_status', 'tumor_only_snv']].copy()

    upset_plot_df = upset_plot_df[upset_plot_df['model_consensus_purity'] > 0.8]
    upset_plot_df = upset_plot_df.drop(columns = 'model_consensus_purity')

    upset_plot_df['loh_shared'] = np.where(upset_plot_df['loh_shared'] > 0.8, 1, 0)
    upset_plot_df['cosine_similarity'] = np.where(upset_plot_df['cosine_similarity'] > 0.9, 1, 0)
    upset_plot_df['wgd_status'] = np.where((upset_plot_df['wgd_status'] == 'Both WGD') | (upset_plot_df['wgd_status'] == 'Both diploid'), 1, 0)
    upset_plot_df['tumor_only_snv'] = np.where(upset_plot_df['tumor_only_snv'] < 0.2, 1, 0)

    matrix_for_class_count = upset_plot_df



    matrix_for_class_count['mismatch_count'] = matrix_for_class_count['loh_shared'] + matrix_for_class_count['cosine_similarity'] + matrix_for_class_count['wgd_status'] + matrix_for_class_count['tumor_only_snv']
    counts = Counter(matrix_for_class_count['mismatch_count'])
    sorted_counts = sorted(counts.items()) 
    fractions = [val / len(matrix_for_class_count) for _, val in sorted_counts]


    barplot_df = pd.DataFrame(sorted_counts)
    barplot_df.columns = ['num_concordant_classes', 'num_samples']

    three_or_more_concordant = barplot_df.loc[barplot_df['num_concordant_classes'] >= 3, 'num_samples'].sum()
    two_concordant = barplot_df.loc[barplot_df['num_concordant_classes'] == 2, 'num_samples'].sum()
    one_concordant = barplot_df.loc[barplot_df['num_concordant_classes'] == 1, 'num_samples'].sum()
    none_concordant = barplot_df.loc[barplot_df['num_concordant_classes'] == 0, 'num_samples'].sum()

    df_for_barplot = pd.DataFrame({
        'three_concordant': [three_or_more_concordant],
        'two_concordant': [two_concordant],
        'one_concordant': [one_concordant],
        'none_concordant': [none_concordant]
    })
    df_for_barplot = df_for_barplot.T
    df_for_barplot.columns = ['num_samples']

    sns.barplot(data = df_for_barplot, x = df_for_barplot.index, y = 'num_samples')
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig("concordance_barplot.pdf", dpi=600, bbox_inches='tight')




#Do the thing
def main():
    sample_list = load_sample_list()
    mutsig_data = load_mutsig_data()
    consensus_data = load_consensus_data()
    snv_consensus_data = load_snv_consensus_data()
    construct_consensus_barplot(sample_list, mutsig_data, consensus_data, snv_consensus_data)




if __name__ == '__main__':
    main()