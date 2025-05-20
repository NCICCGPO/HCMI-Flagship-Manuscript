import matplotlib
import matplotlib.pyplot as plt

import click
import seaborn as sns
import pandas as pd

from gdan_hcmi_tools.plots import palette_wgd

def plot_fig2D(pairwise_info, wgd_palette, fig_path=None):
    plt.rcParams['pdf.fonttype'] = 42
    wgd_columns = ['Both diploid', 'Both WGD', 'Model WGD', 'Tumor WGD']
    with matplotlib.rc_context({'font.family':'Arial'}):
        fig, ax = plt.subplots(figsize=(3.5, 3))
        sns.scatterplot(pairwise_info, x='tumor_consensus_ploidy', y='model_consensus_ploidy', hue='wgd_status', palette=wgd_palette,
                    s=25, alpha=0.8, legend=False)
        handles, labels = ax.get_legend_handles_labels()
        
        #specify order of items in legend
        order = [0, 2, 1, 3]
        
        #add legend to plot
        # ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], 
        #           frameon=False, loc='upper left', bbox_to_anchor=(1,1))
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.plot([2, 2], [0, 10], linestyle='--', lw=1, color='grey', zorder=-1)
        ax.plot([4, 4], [0, 10], linestyle='--', lw=1, color='grey', zorder=-1)
        ax.plot([0, 10], [2, 2], linestyle='--', lw=1, color='grey', zorder=-1)
        ax.plot([0, 10], [4, 4], linestyle='--', lw=1, color='grey', zorder=-1)

        # ax.legend()
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel('Tumor Ploidy')
        ax.set_ylabel('Model Ploidy')
        if fig_path:
            fig.savefig(fig_path, bbox_inches='tight')

@click.command()
@click.option('--pairwise', type=str, help="input pairwise info table")
@click.option('-o', '--output', type=str, help="output figure path for fig2D")
def make_fig2D(pairwise, output):

    pairwise_info = pd.read_csv(pairwise)
    plot_fig2D(pairwise_info, palette_wgd, output)

if __name__ == "__main__":
    make_fig2D()