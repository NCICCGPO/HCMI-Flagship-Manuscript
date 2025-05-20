import click
import pandas as pd


def make_summary_from_notes(notes):
    notes = notes.sort_values('count', ascending=False).reset_index(drop=True)
    summary = (notes[['cancer_type', 'gene', 'count']]
        .groupby(['cancer_type', 'gene'])
        .agg({'count':'sum'})
        .sort_values('count', ascending=False)
        .reset_index())
    return summary

def make_joint_freq_table(summary, counts, hcmi_freq):
    data = []
    counts = counts.set_index('cancer_type')
    for tt in hcmi_freq['cancer_type'].unique():
        tt_stats = summary[summary['cancer_type']==tt].set_index('gene')
        hcmi_stats = hcmi_freq[hcmi_freq['cancer_type']==tt].set_index('gene')
        hcmi_n = int(hcmi_stats['n_samples'].iloc[0])
        for gene, row in hcmi_stats.iterrows():
            if gene == 'TERT': # TERT promoter mutations does not exist in TCGA; unjust comparison
                continue
            hfreq = row['total']
            if tt not in counts.index.values:
                continue
            tcga_n = counts.loc[tt].squeeze()
            tcnt = 0
            if gene in tt_stats.index.values:
                tcnt = tt_stats.loc[gene, 'count']
                tfreq = tcnt / tcga_n
                field = [gene, hfreq, tfreq, tt, hcmi_n, tcga_n]
                data.append(field)
    freqs = pd.DataFrame(data, columns=['gene', 'frequency_HCMI', 'frequency_TCGA', 'cancer_type', 'n_HCMI', 'n_TCGA'])
    return freqs

@click.command()
@click.option('--summary', type=str, help="input TCGA tumor-model comparison summary")
@click.option('--counts', type=str, help="input TCGA tumor type counts from available data")
@click.option('--notes_hcmi', type=str, help="input HCMI tumor-model comparison notes")
@click.option('--cohort_hcmi', type=str, help="input HCMI cohort table")
@click.option('--freq', type=str, help="output frequency comparison data for figure")
def make_fig3A_data(summary, counts, notes_hcmi, cohort_hcmi, freq):
    summary = pd.read_csv(summary)
    counts = pd.read_csv(counts)

    hcmi_obs = pd.read_csv(cohort_hcmi)
    hcmi_obs = hcmi_obs[hcmi_obs['sample_type'].str.startswith('model')]
    hcmi_counts = hcmi_obs['cancer_type'].value_counts()

    hcmi_notes = pd.read_csv(notes_hcmi)
    hcmi_notes = hcmi_notes[hcmi_notes['note']!='tumor-only']
    hcmi_notes.drop('note', axis=1, inplace=True)
    hcmi_summary = make_summary_from_notes(hcmi_notes)
    for rix, row in hcmi_summary.iterrows():
        tt = row['cancer_type'] # tt: tumor type
        _count = row['count']
        _n_models = hcmi_counts.loc[tt]
        _freq = _count / _n_models
        hcmi_summary.loc[rix, 'total'] = _freq
        hcmi_summary.loc[rix, 'n_samples'] = _n_models
    hcmi_freq = hcmi_summary[['gene', 'total', 'cancer_type', 'n_samples']]
    joint_freq = make_joint_freq_table(summary, counts, hcmi_freq)
    joint_freq.to_csv(freq, index=False)

if __name__ == "__main__":
    make_fig3A_data()
