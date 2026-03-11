
# source("OncoMatch/blob-generator.R")

source("../vaxtools/R/interactome_handler.R")
source("../vaxtools/R/utils.R")
# source("../vaxtools/R/cross-species-utils.R")
source("OncoMatch/blob-utils.R")

# create_workspace(run_dir = "blobs-gbm-only-positive-mrs")
create_workspace(run_dir = "blobs-gbm-positive-and-negative-mrs-unbiased")

n_mrs <- 5
run_name <- "gbm"
# run_name <- "brca"
# run_name <- "coad"
# run_name <- "paad"

om_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/oncomatch/gdat-OncoMatch-model2tumor-datafreeze20240311_null-bootstrap_Glioblastoma.rds")
# om_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/oncomatch/gdat-OncoMatch-model2tumor-datafreeze20240311_null-bootstrap_Breast Cancer.rds")
# om_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/oncomatch/gdat-OncoMatch-model2tumor-datafreeze20240311_null-bootstrap_Colorectal Cancer.rds")
# om_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/oncomatch/gdat-OncoMatch-model2tumor-datafreeze20240311_null-bootstrap_Pancreatic Cancer.rds")
om_data <- as_tibble(om_data)

om_data <- om_data %>% arrange(pval)

selected_model_id_good_match <- om_data$model.id[1]
selected_tumor_id_good_match <- om_data$tumor.id[1]
selected_model_id_bad_match <- om_data$model.id[nrow(om_data)]
selected_tumor_id_bad_match <- om_data$tumor.id[nrow(om_data)]

# ## GBM Cases for Model drifts (Mesenchymal vs Proneural)
# selected_model_id_good_match <- "HCM-BROD-0415-C71-85R"
# selected_tumor_id_good_match <- "HCM-BROD-0415-C71-02B"
# selected_model_id_bad_match <- "HCM-BROD-0002-C71-85A"
# selected_tumor_id_bad_match <- "HCM-BROD-0002-C71-01A"

# ## Some COAD Cases Heeju suggested
# selected_model_id_good_match <- "HCM-STAN-1111-C19-85B"
# selected_tumor_id_good_match <- "HCM-STAN-1111-C19-01A"
# # selected_model_id_bad_match <- "HCM-CSHL-0255-C18-85F"
# # selected_tumor_id_bad_match <- "HCM-CSHL-0255-C18-01A"
# selected_model_id_bad_match <- "HCM-CSHL-0258-C18-85D"
# selected_tumor_id_bad_match <- "HCM-CSHL-0258-C18-01A"

## Protein Activity Signature Selection ----
pas_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/pmat/pmat-HCMI-dataFreeze2-20240311_center-hcmi-GBM-LGG-noIntegration-tumor-model_regulGBM_aracne3-mi00-size100_narnea.rds")
# pas_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/pmat/pmat-HCMI-dataFreeze2-20240311_center-hcmi-BRCA-noIntegration-tumor-model_regulBRCA_aracne3-mi00-size100_narnea.rds")
# pas_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/pmat/pmat-HCMI-dataFreeze2-20240311_center-hcmi-COAD-READ-noIntegration-tumor-model_regulCOAD_aracne3-mi00-size100_narnea.rds")
# pas_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/pmat/pmat-HCMI-dataFreeze2-20240311_center-hcmi-PAAD-noIntegration-tumor-model_regulPAAD_aracne3-mi00-size100_narnea.rds")
dim(pas_data)

pas_signature_model_good_match <- pas_data[,selected_model_id_good_match]
pas_signature_tumor_good_match <- pas_data[,selected_tumor_id_good_match]
pas_signature_model_bad_match <- pas_data[,selected_model_id_bad_match]
pas_signature_tumor_bad_match <- pas_data[,selected_tumor_id_bad_match]

# tail(sort(pas_signature_model_good_match))
# pas_signature_tumor_good_match[ names(tail(sort(pas_signature_model_good_match))) ]

# selected_model_mrs_good_match <- getTopMRs( pas_data[,selected_model_id_good_match] , n_top = n_mrs , onlyPositive = T)
# selected_tumor_mrs_good_match <- getTopMRs( pas_data[,selected_tumor_id_good_match] , n_top = n_mrs , onlyPositive = T)
# selected_model_mrs_bad_match <- getTopMRs( pas_data[,selected_model_id_bad_match] , n_top = n_mrs , onlyPositive = T)
# selected_tumor_mrs_bad_match <- getTopMRs( pas_data[,selected_tumor_id_bad_match] , n_top = n_mrs , onlyPositive = T)
selected_model_mrs_good_match <- getTopMRs( pas_data[,selected_model_id_good_match] , n_top = n_mrs , onlyPositive = F)
selected_tumor_mrs_good_match <- getTopMRs( pas_data[,selected_tumor_id_good_match] , n_top = n_mrs , onlyPositive = F)
selected_model_mrs_bad_match <- getTopMRs( pas_data[,selected_model_id_bad_match] , n_top = n_mrs , onlyPositive = F)
selected_tumor_mrs_bad_match <- getTopMRs( pas_data[,selected_tumor_id_bad_match] , n_top = n_mrs , onlyPositive = F)

all_mrs_good_match <- c(selected_model_mrs_good_match,selected_tumor_mrs_good_match)
all_mrs_bad_match <- c(selected_model_mrs_bad_match,selected_tumor_mrs_bad_match)
# all_mrs_good_match <- c("STAT3","CEBPB","CEBPD","CD44","S100A16","WWTR1","TNFRSF1A","ANXA1","ANXA2",
# 						"SOX2","MYCN","OLIG2","DLL3","PDGFRA","STMN4", "DLL1", "ASCL1", "DLX5")
# all_mrs_bad_match <- all_mrs_good_match


## Gene Expression Signature Selection ----
ges_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/signatures/gesig-hcmi-centerByTissue_GBM.rds")
# ges_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/signatures/gesig-hcmi-centerByTissue_BRCA.txt")
# ges_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/signatures/gesig-hcmi-centerByTissue_COAD.rds")
# ges_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/signatures/gesig-hcmi-centerByTissue_PAAD.rds")
dim(ges_data)

ges_signature_model_good_match <- ges_data[,selected_model_id_good_match]
ges_signature_tumor_good_match <- ges_data[,selected_tumor_id_good_match]
ges_signature_model_bad_match <- ges_data[,selected_model_id_bad_match]
ges_signature_tumor_bad_match <- ges_data[,selected_tumor_id_bad_match]

## Network manipulation: Generated edges and nodes table ----
net_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/networks/regulon-TCGA-GBM_aracne3-mi0-regSize100.rds")
# net_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/networks/regulon-TCGA-COAD_aracne3-mi0-regSize100.rds")
# net_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/networks/regulon-TCGA-PAAD_aracne3-mi0-regSize100.rds")
net_data <- net_data$GBM
# net_data <- net_data$PAAD
# net_data <- net_data$COAD

##########################################
## Plotting Good Match MRs over Blobs ----
##########################################
## MODEL --
res <- generateNodesAndEdges(net=net_data,mrs_list=all_mrs_good_match)
res <- integrateExpressionAndActivity(edges_table = res$edges , nodes_table = res$nodes ,
									ges_signature_vector = ges_signature_model_good_match , 
									pas_signature_vector = pas_signature_model_good_match)
nodes_model_good_match <- res$nodes
nodes_model_good_match[ nodes_model_good_match$type == "TF" , ]
filename <- file.path(reports.dir, paste0("blobs-",run_name,"-good-match",
										  "-model-",selected_model_id_good_match,"-tumor-",selected_tumor_id_good_match,
										  "-modelSignature.pdf"))
plotBlobs( edges_table = res$edges , nodes_table = res$nodes , blob_filename = filename)

## TUMOR --
res <- generateNodesAndEdges(net=net_data,mrs_list=all_mrs_good_match)
res <- integrateExpressionAndActivity(edges_table = res$edges , nodes_table = res$nodes ,
									  ges_signature_vector = ges_signature_tumor_good_match , 
									  pas_signature_vector = pas_signature_tumor_good_match)
nodes_tumor_good_match <- res$nodes
nodes_tumor_good_match[ nodes_tumor_good_match$type == "TF" , ]
filename <- file.path(reports.dir, paste0("blobs-",run_name,"-good-match",
										  "-model-",selected_model_id_good_match,"-tumor-",selected_tumor_id_good_match,
										  "-tumorSignature.pdf"))
plotBlobs( edges_table = res$edges , nodes_table = res$nodes , blob_filename = filename)

#########################################
## Plotting Bad Match MRs over Blobs ----
#########################################
## MODEL --
res <- generateNodesAndEdges(net=net_data,mrs_list=all_mrs_bad_match)
res <- integrateExpressionAndActivity(edges_table = res$edges , nodes_table = res$nodes ,
									  ges_signature_vector = ges_signature_model_bad_match , 
									  pas_signature_vector = pas_signature_model_bad_match)
nodes_model_bad_match <- res$nodes
nodes_model_bad_match[ nodes_model_bad_match$type == "TF" , ]
filename <- file.path(reports.dir, paste0("blobs-",run_name,"-bad-match",
										  "-model-",selected_model_id_bad_match,"-tumor-",selected_tumor_id_bad_match,
										  "-modelSignature.pdf"))
plotBlobs( edges_table = res$edges , nodes_table = res$nodes , blob_filename = filename)

## TUMOR --
res <- generateNodesAndEdges(net=net_data,mrs_list=all_mrs_bad_match)
res <- integrateExpressionAndActivity(edges_table = res$edges , nodes_table = res$nodes ,
									  ges_signature_vector = ges_signature_tumor_bad_match , 
									  pas_signature_vector = pas_signature_tumor_bad_match)
nodes_tumor_bad_match <- res$nodes
nodes_tumor_bad_match[ nodes_tumor_bad_match$type == "TF" , ]
filename <- file.path(reports.dir, paste0("blobs-",run_name,"-bad-match",
										  "-model-",selected_model_id_bad_match,"-tumor-",selected_tumor_id_bad_match,
										  "-tumorSignature.pdf"))
plotBlobs( edges_table = res$edges , nodes_table = res$nodes , blob_filename = filename)






