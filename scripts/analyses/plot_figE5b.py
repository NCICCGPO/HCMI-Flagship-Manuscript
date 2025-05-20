import matplotlib
import matplotlib.pyplot as plt

import click
import numpy as np
import pandas as pd

@click.command()
@click.option('--freq', type=str, help="input frequency comparison data for figure")
@click.option('--min_samples_per_cancer_type', type=int, help="minimum number of samples per cancer type to be included in the plot", default=10, show_default=True)
@click.option('--n_top_genes_to_show', type=int, help="number of most prevalent genes to show per cancer type (0 shows all genes)", default=0, show_default=True)
@click.option('--figwidth', type=int, help="figure width in inches", default=16, show_default=True)
@click.option('--figheight', type=int, help="figure height in inches", default=4, show_default=True)
@click.option('--nrow', type=int, help="number of rows for the figure", default=2, show_default=True)
@click.option('--fig_path', type=str, help="output fig path")
def plot_fig3A(freq, min_samples_per_cancer_type, n_top_genes_to_show, figwidth, figheight, nrow, fig_path):
    freqs = pd.read_csv(freq)
    freqs.sort_values(['n_HCMI', 'frequency_HCMI', 'gene'], ascending=[False, False, True], inplace=True)
    freqs = freqs[(freqs['n_HCMI'] >= min_samples_per_cancer_type) & (freqs['n_TCGA'] >= min_samples_per_cancer_type)]
    tts = freqs['cancer_type'].unique().tolist() # tumor types
    cohort_colors = {'HCMI':'#149bd5', 'TCGA':'lightgrey'}
    width = 0.4
    n_top_genes_to_show = 5
    ncol = len(tts) // nrow + len(tts) % nrow
    fig_y = 1.05
    if figheight > 5:
        fig_y = 1.0

    with matplotlib.rc_context({'font.family':'Arial'}):
        fig, axes = plt.subplots(nrow, ncol, figsize=(figwidth, figheight), gridspec_kw={'hspace':1.3, 'wspace':1})
        fig.suptitle('Driver events', y=fig_y)
        for ix, tt in enumerate(tts):
            tt_freq = freqs[freqs['cancer_type']==tt].iloc[:n_top_genes_to_show]
            ax = axes.flatten()[ix]
            n_hcmi = int(tt_freq['n_HCMI'].iloc[0])
            n_tcga = int(tt_freq['n_TCGA'].iloc[0])
            for cix, cohort in enumerate(['HCMI', 'TCGA']):
                color = cohort_colors[cohort]
                col = f'frequency_{cohort}'
                xs = np.arange(tt_freq.shape[0])
                heights = tt_freq[col]
                ax.bar(x=xs + cix*width, height=heights, width=width, color=color, label=cohort)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            ax.text(x=(xmin+xmax)/2 * 0.8, y=1.15 * ymax, s=f'(HCMI: {n_hcmi:d}, TCGA: {n_tcga:d})', ha='center', va='center')
            ax.set_title(f'{tt}\n', fontsize=12)
            ax.set_xticks(np.arange(tt_freq.shape[0]) + width/2)
            ax.set_xticklabels(tt_freq['gene'], rotation=90)
            ax.spines[['right', 'top']].set_visible(False)
            show_legend = ix == len(tts)-1
            if show_legend:
                ax.legend(loc='upper left', bbox_to_anchor=(1,1), frameon=False)
        for jx in range(ix+1, nrow*ncol):
            ax = axes.flatten()[jx]
            ax.axis('off')
    fig.savefig(fig_path, bbox_inches='tight')

if __name__ == "__main__":
    plot_fig3A()
