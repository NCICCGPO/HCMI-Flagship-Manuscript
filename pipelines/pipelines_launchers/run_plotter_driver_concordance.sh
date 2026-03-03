python scripts/analyses/plot_fig2C.py \
    --pairwise data/driver-concordance/tumor_model_wgd_features.csv \
    -o1 plots/fig2C_1.pdf \
    -o2 plots/fig2C_2.pdf
echo "[LOG] fig2C done"
python scripts/analyses/plot_fig2D.py \
    --pairwise data/driver-concordance/tumor_model_wgd_features.csv \
    -o plots/fig2D.pdf
echo "[LOG] fig2D done"
python scripts/analyses/plot_fig2E.py \
    --summary data/driver-concordance/concordance_summary.csv \
    --counts data/driver-concordance/tumor_sample_counts.csv \
    --cn_drivers_path data/driver-concordance/cnv_drivers.yaml \
    --fig_path plots/fig2E.pdf
echo "[LOG] fig2E done"
python scripts/analyses/make_data_for_figE5b.py \
    --summary data/driver-concordance/tcga_concordance_summary.csv \
    --counts data/driver-concordance/tcga_tumor_sample_counts.csv \
    --notes_hcmi data/driver-concordance/concordance_per_mutation_type.csv \
    --cohort_hcmi data/driver-concordance/snv_indel_cnv_pair_information.csv \
    --freq data/driver-concordance/event_frequency_hcmi_tcga.csv
python scripts/analyses/plot_figE5b.py \
    --freq data/driver-concordance/event_frequency_hcmi_tcga.csv \
    --fig_path plots/extended_fig5b.pdf \
    --min_samples_per_cancer_type 1 \
    --n_top_genes_to_show 10
echo "[LOG] figE5b done"
