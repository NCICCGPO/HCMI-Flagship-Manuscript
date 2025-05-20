import click
import yaml
import pandas as pd
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from gdan_hcmi_tools.plots import palette_concordance


def plot_concordance(concordances, palette, title='', figsize=(15, 3), top_n=50, 
        ax=None, show_legend=False, show_xlabel=True, genes=None, return_fig=False,):
    plt.rcParams['pdf.fonttype'] = 42
    plot_data = concordances.iloc[:top_n].copy()
    genes2add = []
    if genes:
        for gene, row in concordances.iterrows():
            if gene in genes and gene not in plot_data.index:
                genes2add.append(gene)
    plot_data = pd.concat([plot_data, concordances.loc[genes2add]])
    plot_data.drop('total', axis=1, inplace=True)
    plot_data = plot_data[['both', 'model-only', 'tumor-only']]
    plot_data.rename(columns={'both':'Intersecting', 'model-only':'Model only', 'tumor-only':'Tumor only'}, inplace=True)
    color_order = ['Intersecting', 'Model only', 'Tumor only']
    with matplotlib.rc_context({'font.family':'Arial'}):
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=figsize)
        else:
            fig = ax.get_figure()
        plot_data.plot(kind='bar', stacked=True, ax=ax, color=palette)
        if genes:
            new_xticklabels = []
            for ix, xticklabel in enumerate(ax.get_xticklabels()):
                text = xticklabel.get_text()
                color = 'black'
                if text in genes and ix >= 5:
                    color = 'red'
                xticklabel.set_c(color)
                new_xticklabels.append(xticklabel)
            ax.set_xticklabels(new_xticklabels)
        if title:
            ax.set_title(title)
        if show_legend:
            handles = [mpatches.Patch(label=x, color=palette[x]) for x in color_order]
            ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
        else:
            ax.get_legend().remove()
        ax.spines[['right', 'top']].set_visible(False)
        if not show_xlabel:
            ax.set_xlabel('')
    if return_fig:
        return fig


def plot_fig2E(summary, counts, cn_drivers, fig_path, palette, min_count=5):
    nrow, ncol = 3, 6
    counts = counts[counts['count'] >= min_count].copy()
    tts = counts.index.tolist() # tumor types
    plt.rcParams['pdf.fonttype'] = 42
    with matplotlib.rc_context({'font.family':'Arial'}):
        fig, axes = plt.subplots(nrow, ncol, figsize=(12, 7), gridspec_kw={'hspace':1.2, 'wspace':.5})
        fig.suptitle('Driver events')
        for ix in range(nrow * ncol):
            if ix < counts.shape[0]:
                tt = tts[ix] # tumor type
                tt_count = counts.loc[tt, 'count']
                show_legend = (tt=='CHOL')
                tt_df = (summary[summary['cancer_type']==tt]
                    .drop('cancer_type', axis=1)
                    .set_index('gene'))
                tt_df = tt_df.copy() / tt_count # normalize to frequency
                ax = axes.flatten()[ix]
                if tt_df.shape[0] > 0:
                    title = f'{tt} (n={tt_count})'
                    genes = None
                    if tt in cn_drivers:
                        genes = cn_drivers[tt]
                    plot_concordance(tt_df, palette, title=title, figsize=(3, 1), top_n=5, ax=ax, 
                        show_legend=show_legend, show_xlabel=False, genes=genes)
            else:
                ax = axes.flatten()[ix]
                ax.axis('off')
    fig.savefig(fig_path, bbox_inches='tight')


@click.command()
@click.option('--summary', type=str, help="input tumor-model comparison summary")
@click.option('--counts', type=str, help="input tumor type counts from available data")
@click.option('--cn_drivers_path', type=str, help="input cancer type specific cnv drivers")
@click.option('--fig_path', type=str, help="output figure path")
def make_figure2E(summary, counts, cn_drivers_path, fig_path):
    pivot_summary = pd.read_csv(summary)
    tumor_type_counts = pd.read_csv(counts, index_col=0)

    cn_drivers = yaml.safe_load(open(cn_drivers_path))
    plot_fig2E(pivot_summary, tumor_type_counts, cn_drivers, fig_path, palette_concordance)

if __name__ == "__main__":
    make_figure2E()
