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

sample_list = load_sample_list()


#Load in the signature matrix
def load_signature_matrix():
    signatures = pd.read_hdf('nmf_output.h5', 'H')
    signatures.index = signatures.index.str.split('-', n = 5).str[:5].str.join('-')
    signatures.drop(columns=['max', 'max_id', 'max_norm'], inplace=True)

    row_sums = signatures.sum(axis = 1)
    signatures_normalized = signatures.div(row_sums, axis = 0)
    return signatures_normalized



#Define a function to compute cosine similarity
def compute_cosine_similarity(sample_list, signatures_normalized):
    random.seed(42)
    sample_list = sample_list.copy()
    signatures_normalized = signatures_normalized.copy()

    #sometimes there are samples in the pair list that we don't have signature scores for. Let's just go ahead and remove those pairs
    sample_list = sample_list[
        (sample_list['model'].isin(signatures_normalized.index)) | (sample_list['tumor'].isin(signatures_normalized.index))
    ]

    models = list(sample_list['model'])
    tumors = list(sample_list['tumor'])

    #First get a tuples of tumor-model pairs
    tumor_model_pair_tuple = list(zip(models, tumors))
    random_pairs = list(zip(models, random.sample(tumors, len(tumors))))

    
    #Get lists of cosine similarity scores
    def get_cosine_similarity_scores(pairs):
        cosine_similarity_scores = []
        for pair in pairs:

            if (pair[0] in signatures_normalized.index) and (pair[1] in signatures_normalized.index):
                model_scores = signatures_normalized.loc[pair[0]]
                tumor_scores = signatures_normalized.loc[pair[1]]
                cosine_similarity = compute_tumor_model_cosine_similarity(model_scores, tumor_scores)
                cosine_similarity_scores.append(cosine_similarity)
            else:
                continue
        
        return cosine_similarity_scores
    


    #now do the same for random pairs within tumor type
    def compute_tumor_type_restricted(models, tumors, pair_cosine_similarity_list, signatures_normalized, comparison):
        tested = []
        cosine_similarity_scores = []

        while len(tested) < len(pair_cosine_similarity_list):
            random_model = random.choice(models)
            random_tumor = random.choice(tumors)
            model_type = random_model.split("-")[3]
            tumor_type = random_tumor.split("-")[3]

            if (random_model, random_tumor) in tested:
                continue

            if model_type != tumor_type and comparison == 'same':
                continue

            if model_type == tumor_type and comparison == 'different':
                continue

            if (random_model in signatures_normalized.index) and (random_tumor in signatures_normalized.index):
                model_scores = signatures_normalized.loc[random_model]
                tumor_scores = signatures_normalized.loc[random_tumor]
                cosine_similarity = compute_tumor_model_cosine_similarity(model_scores, tumor_scores)
                tested.append((random_model, random_tumor))
                cosine_similarity_scores.append(cosine_similarity)

        return cosine_similarity_scores


    #now run it all
    pair_cosine_similarity_list = get_cosine_similarity_scores(tumor_model_pair_tuple)
    random_cosine_similarity_list = get_cosine_similarity_scores(random_pairs)
    within_tumor_type_cosine_similarity = compute_tumor_type_restricted(models, tumors, pair_cosine_similarity_list, signatures_normalized, comparison = 'same')
    not_within_tumor_type_cosine_similarity = compute_tumor_type_restricted(models, tumors, pair_cosine_similarity_list, signatures_normalized, comparison = 'different')
    
    #Construct a final output matrix
    score_df = pd.DataFrame({
        'pair': pair_cosine_similarity_list,
        'random': random_cosine_similarity_list,
        'random_same_tumor_type': within_tumor_type_cosine_similarity,
        'random_different_tumor_type': not_within_tumor_type_cosine_similarity
    })

    return score_df







#Define a function to construct the violin plot
def constuct_fig2_mutation_sig_violin(score_df):

    score_df_melted = score_df.melt()


    #Plot the violin plot
    sns.violinplot(
        x = 'variable', 
        y = 'value', 
        data = score_df_melted, 
        inner = None,
        color = 'grey',
        linewidth = 0.5,
        alpha = 0.5,
        width = 0.8
    )

    #Plot the strip plot
    sns.stripplot(x = 'variable', 
                y = 'value', 
                data = score_df_melted, 
                alpha = 1, 
                edgecolor = 'darkred',
                size = 2,
                color = "darkred", 
                linewidth = 0.5)
    plt.axhline(y=0.9, color='black', linestyle='dotted', linewidth=0.5)

    #Change some figure parameters
    plt.xticks(rotation=45)

    #Export the figure
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig("signature_overlapping_violin_with_controls.pdf", dpi=600, bbox_inches='tight')

constuct_fig2_mutation_sig_violin(score_df)



#Do the thing
def main():
    sample_list = load_sample_list()
    signatures_normalized = load_signature_matrix()
    cosine_similarity_scores = compute_cosine_similarity(sample_list, signatures_normalized)
    score_df = compute_cosine_similarity(sample_list, signatures_normalized)
    constuct_fig2_mutation_sig_violin(score_df)


if __name__ == '__main__':
    main()