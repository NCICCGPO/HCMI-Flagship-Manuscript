## Custom functions and colors for DNA methylation analysis
## Date: 2025-05-30

# Fig3a, Fig4c, and Fig4e
probe_stats <- function(data = dt, me_cutoff = 0.3) {
  probe_id <- data[, ProbeID]
  num_all <- ncol(data) - 1
  num_na <- rowSums(is.na(data[, .SD, .SDcols = -"ProbeID"]))
  num_data <- rowSums(!is.na(data[, .SD, .SDcols = -"ProbeID"]))
  num_meth <- rowSums(data[, .SD, .SDcols = -"ProbeID"] >= me_cutoff,
                      na.rm = TRUE)
  probe_sds <- data[, rowSds(as.matrix(.SD), na.rm = TRUE),
                    .SDcols = -"ProbeID"]
  perc_na <- (num_na / num_all) * 100
  perc_meth <- (num_meth / num_data) * 100
  return(data.table(probe_id,
                    num_all,
                    num_na,
                    num_data,
                    num_meth,
                    perc_na,
                    perc_meth,
                    probe_sds))
}

## Fig4c and Extended Data Fig.5c
theme_fig4c <- theme_bw() +
  theme(
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1)

## Fig4c and Extended Data Fig.5c
hcmi_tcga_colors <-
  c(
    "ACC" = "#C1A72F",
    "AITL" = "#FF938B",
    "AMPV" = "#E8C51D",
    "BLCA" = "#FAD2D9",
    "BRCA" = "#ED2891",
    "CESC" = "#F6B667",
    "CHOL" = "#104A7F",
    "COAD" = "#9EDDF9",
    "COTA" = "#F89420",
    "CTCL" = "#7E1918",
    "DCIS" = "#ff4d75",
    "DEST" = "#FF6347",
    "DLBC" = "#3953A4",
    "ESCA" = "#007EB5",
    "GBAD" = "#8A2BE2",
    "GBM" = "#B2509E",
    "HNSC" = "#97D1A9",
    "IPMN" = "#542C88",
    "KICH" = "#ED1C24",
    "KIRC" = "#F8AFB3",
    "KIRP" = "#EA7075",
    "LAML" = "#754C29",
    "LGG" = "#D49DC7",
    "LIHC" = "#CACCDB",
    "LUAD" = "#D3C3E0",
    "LUSC" = "#A084BD",
    "MB" = "#3953A4",
    "MESO" = "#542C88",
    "MISC" = "#00CED1",
    "NB" = "#20B2AA",
    "NHL" = "#BE1E2D",
    "OV" = "#D97D25",
    "PAAD" = "#6E7BA2",
    "PCPG" = "#E8C51D",
    "PRAD" = "#7E1918",
    "PTCL" = "#BC8279",
    "READ" = "#DAF1FC",
    "RMC" = "#FF4500",
    "SARC" = "#00A99D",
    "SCLC" = "#8B4513",
    "SESP" = "#D49DC7",
    "SIAC" = "#FFD700",
    "SIAD" = "#009444",
    "SKCM" = "#BBD642",
    "SSCC" = "#4682B4",
    "STAD" = "#00AEEF",
    "TGCT" = "#BE1E2D",
    "THCA" = "#F9ED32",
    "THYM" = "#CEAC8F",
    "TPLL" = "#EA7075",
    "TVA" = "#C1A72F",
    "UCEC" = "#FBE3C7",
    "UCS" = "#F89420",
    "UVM" = "#BBD642",
    "WT" = "#754C29"
  )

color_map_border <- c("TCGA_TUMOR" = NA,
                      "TARGET_TUMOR" = NA,
                      "HCMI_TUMOR" = "black",
                      "HCMI_MODEL" = "black")

shape_map <- c("TCGA_TUMOR" = NA,
               "TARGET_TUMOR" = NA,
               "HCMI_TUMOR" = 24,
               "HCMI_MODEL" = 21)

alpha_map <- c("TCGA_TUMOR" = 0.2,
               "TARGET_TUMOR" = 0.2,
               "HCMI_TUMOR" = 1,
               "HCMI_MODEL" = 1)


## Fig4e
project_colors <- c(
  "TCGA" = "#CCCCCC",
  "HCMI" = "#000000"
)

sample_type_colors <-
  c("Tumor" = "#FFFFFF",
    "Model_3D" = "#0000FF",
    "Model_2D" = "#FF0000")

tmp_subtype_colors <-
  c("COADREAD_1" = "#E6194B",
    "COADREAD_2" = "#3CB44B",
    "COADREAD_3" = "#FFE119",
    "COADREAD_4" = "#4363D8",
    "COADREAD_5" = "#F58231",
    "COADREAD_6" = "#911EB4",
    "GEA_1" = "#E6194B",
    "GEA_2" = "#3CB44B",
    "GEA_3" = "#FFE119",
    "GEA_4" = "#4363D8",
    "GEA_5" = "#F58231",
    "GEA_6" = "#911EB4",
    "GEA_7" = "#42D4F4",
    "LGGGBM_2" = "#3CB44B",
    "LGGGBM_3" = "#FFE119",
    "LGGGBM_4" = "#4363D8",
    "LGGGBM_5" = "#F58231",
    "LGGGBM_6" = "#911EB4",
    "LGGGBM_7" = "#42D4F4",
    "LGGGBM_1" = "#E6194B",
    "LGGGBM_2" = "#3CB44B",
    "LGGGBM_3" = "#FFE119",
    "LGGGBM_4" = "#4363D8",
    "LGGGBM_5" = "#F58231",
    "LGGGBM_6" = "#911EB4",
    "LGGGBM_7" = "#42D4F4",
    "PAAD_1" = "#E6194B",
    "PAAD_2" = "#3CB44B")

col_fun_turbo100 <- circlize::colorRamp2(seq(from = 0, to = 1,
                                             length.out = 100),
                                         viridis::turbo(100))
