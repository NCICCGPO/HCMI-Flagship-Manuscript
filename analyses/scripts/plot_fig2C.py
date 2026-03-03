import matplotlib
import matplotlib.pyplot as plt

import click
import seaborn as sns
import pandas as pd

from gdan_hcmi_tools.plots import palette_wgd

def plot_fig2C_per_cancer_type(pairwise_info, wgd_palette, fig_path):
    plt.rcParams['pdf.fonttype'] = 42
    wgd_columns = ['Both diploid', 'Both WGD', 'Model WGD', 'Tumor WGD']
    with matplotlib.rc_context({'font.family':'Arial'}):
        plt.figure(figsize=(2, 4))
        ax = plt.gca()
        wgd_counts = pairwise_info.groupby(['cancer_type', 'wgd_status']).size().unstack(fill_value=0)
        wgd_sum = wgd_counts.sum(axis=1).sort_values()
        others_cts = wgd_sum[wgd_sum < 5].index.tolist()
        others_cnt = wgd_sum.loc[others_cts].sum()
        wgd_sum = wgd_sum[wgd_sum>5]
        high_cts = wgd_sum.index.tolist()
        others_df = pd.DataFrame(wgd_counts.loc[others_cts].sum()).T
        others_df.index = ['Others']
        wgd_counts = pd.concat([wgd_counts.loc[high_cts], others_df])
        wgd_sum.loc['Others'] = others_cnt
        plot_order = wgd_sum.index.tolist()
        plot_order = ['Others'] + [s for s in plot_order if s != 'Others']
        wgd_counts = wgd_counts.loc[plot_order]
        wgd_sum = wgd_sum.loc[plot_order]
        yticklabels = []
        for ct, cnt in wgd_sum.items():
            yticklabel = f'{ct} (n={cnt})'
            yticklabels.append(yticklabel)
        wgd_matrix = (wgd_counts.T / wgd_counts.sum(axis=1)).T[wgd_columns]
        wgd_matrix.plot.barh(ax=ax, stacked=True, color=wgd_palette)
        ax.set_yticklabels(yticklabels)
        sns.despine(ax=ax)
        ax.set_xlabel('Frac. pairs\nWGD status')
        sns.move_legend(ax, loc=(0., 1.03), title='WGD Status', fontsize=8, title_fontsize=10, frameon=False)
        ax.set_ylabel('Cancer type', rotation=90, ha='center', va='center')
        fig = ax.get_figure()
        fig.savefig(fig_path, bbox_inches='tight')

def plot_fig2C_aggregated(pairwise_info, wgd_palette, fig_path):
    plt.rcParams['pdf.fonttype'] = 42
    col_order = ['Both diploid', 'Both WGD', 'Model WGD', 'Tumor WGD']
    with matplotlib.rc_context({'font.family':'Arial'}):
        fig, ax = plt.subplots(figsize=(1, 4))
        wgd_matrix = pairwise_info.assign(placeholder=1).groupby(['placeholder', 'wgd_status']).size().unstack(fill_value=0)
        wgd_matrix = (wgd_matrix.T / wgd_matrix.sum(axis=1)).T
        wgd_matrix.columns = col_order
        wgd_matrix.plot.bar(ax=ax, stacked=True, color=wgd_palette)
        sns.despine(ax=ax)
        ax.set_ylabel('Frac. pairs\nWGD status')
        ax.set_xlabel('')
        ax.set_xticks([])
        sns.move_legend(ax, loc=(0., 1.03), title='WGD Status', fontsize=8, title_fontsize=10, frameon=False)
        fig.savefig(fig_path, bbox_inches='tight')

@click.command()
@click.option('--pairwise', type=str, help="input pairwise info table")
@click.option('-o1', '--fig_path_1', type=str, help="output figure path for fig2C-1 [per cancer type]")
@click.option('-o2', '--fig_path_2', type=str, help="output figure path for fig2C-2 [aggregated]")
def plot_fig2C(pairwise, fig_path_1, fig_path_2):

    pairwise_info = pd.read_csv(pairwise)
    plot_fig2C_per_cancer_type(pairwise_info, palette_wgd, fig_path_1)
    plot_fig2C_aggregated(pairwise_info, palette_wgd, fig_path_2)

if __name__ == "__main__":
    plot_fig2C()