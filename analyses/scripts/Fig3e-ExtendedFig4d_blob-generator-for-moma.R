
# source("OncoMatch/blob-generator-for-moma.R")

source("../vaxtools/R/interactome_handler.R")
source("../vaxtools/R/utils.R")
# source("../vaxtools/R/cross-species-utils.R")
source("OncoMatch/blob-utils.R")

create_workspace(run_dir = "blobs-moma-gbm")

n_mrs <- 5
run_name <- "GBM"

## Loading MOMA data ----
moma_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/moma/tcga-moma-5cohorts-top10MR_sig-centerByTissue.rds")
moma_data <- tibble2matrix(moma_data)
index <- grepl( run_name , rownames(moma_data) )
all_mrs <- moma_data[index,] %>% as.vector()

om_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/oncomatch/gdat-OncoMatch-model2tumor-datafreeze20240311_null-bootstrap_Glioblastoma.rds")
om_data <- as_tibble(om_data)

om_data <- om_data %>% arrange(pval)

selected_model_id_good_match <- om_data$model.id[1]
selected_tumor_id_good_match <- om_data$tumor.id[1]
selected_model_id_bad_match <- om_data$model.id[nrow(om_data)]
selected_tumor_id_bad_match <- om_data$tumor.id[nrow(om_data)]

## Protein Activity Signature Selection ----
pas_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/pmat/pmat-HCMI-dataFreeze2-20240311_center-hcmi-GBM-LGG-noIntegration-tumor-model_regulGBM_aracne3-mi00-size100_narnea.rds")
dim(pas_data)

pas_signature_model_good_match <- pas_data[,selected_model_id_good_match]
pas_signature_tumor_good_match <- pas_data[,selected_tumor_id_good_match]
pas_signature_model_bad_match <- pas_data[,selected_model_id_bad_match]
pas_signature_tumor_bad_match <- pas_data[,selected_tumor_id_bad_match]

# tail(sort(pas_signature_model_good_match))
# pas_signature_tumor_good_match[ names(tail(sort(pas_signature_model_good_match))) ]

selected_model_mrs_good_match <- getTopMRs( pas_data[,selected_model_id_good_match] , n_top = n_mrs , onlyPositive = T)
selected_tumor_mrs_good_match <- getTopMRs( pas_data[,selected_tumor_id_good_match] , n_top = n_mrs , onlyPositive = T)
selected_model_mrs_bad_match <- getTopMRs( pas_data[,selected_model_id_bad_match] , n_top = n_mrs , onlyPositive = T)
selected_tumor_mrs_bad_match <- getTopMRs( pas_data[,selected_tumor_id_bad_match] , n_top = n_mrs , onlyPositive = T)

all_mrs_good_match <- c(selected_model_mrs_good_match,selected_tumor_mrs_good_match)
all_mrs_bad_match <- c(selected_model_mrs_bad_match,selected_tumor_mrs_bad_match)

## Gene Expression Signature Selection ----
ges_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/signatures/gesig-hcmi-centerByTissue_GBM.rds")
dim(ges_data)

ges_signature_model_good_match <- ges_data[,selected_model_id_good_match]
ges_signature_tumor_good_match <- ges_data[,selected_tumor_id_good_match]
ges_signature_model_bad_match <- ges_data[,selected_model_id_bad_match]
ges_signature_tumor_bad_match <- ges_data[,selected_tumor_id_bad_match]

## Network manipulation: Generated edges and nodes table ----
net_data <- readRDS("~/Clouds/Dropbox/Data/hcmi/processed/HCMI_shared with Alessandro/networks/regulon-TCGA-GBM_aracne3-mi0-regSize100.rds")
net_data <- net_data$GBM

##########################################
## Plotting Good Match MRs over Blobs ----
##########################################
## MODEL --
res <- generateNodesAndEdges(net=net_data,mrs_list=all_mrs)
res <- integrateExpressionAndActivity(edges_table = res$edges , nodes_table = res$nodes ,
									  ges_signature_vector = ges_signature_model_good_match , 
									  pas_signature_vector = pas_signature_model_good_match)
nodes_model_good_match <- res$nodes
nodes_model_good_match[ nodes_model_good_match$type == "TF" , ]
filename <- file.path(reports.dir, paste0("blobs-",run_name,"-model-good-match.pdf"))
plotBlobs( edges_table = res$edges , nodes_table = res$nodes , blob_filename = filename)

## TUMOR --
res <- generateNodesAndEdges(net=net_data,mrs_list=all_mrs)
res <- integrateExpressionAndActivity(edges_table = res$edges , nodes_table = res$nodes ,
									  ges_signature_vector = ges_signature_tumor_good_match , 
									  pas_signature_vector = pas_signature_tumor_good_match)
nodes_tumor_good_match <- res$nodes
nodes_tumor_good_match[ nodes_tumor_good_match$type == "TF" , ]
filename <- file.path(reports.dir, paste0("blobs-",run_name,"-tumor-good-match.pdf"))
plotBlobs( edges_table = res$edges , nodes_table = res$nodes , blob_filename = filename)

#########################################
## Plotting Bad Match MRs over Blobs ----
#########################################
## MODEL --
res <- generateNodesAndEdges(net=net_data,mrs_list=all_mrs)
res <- integrateExpressionAndActivity(edges_table = res$edges , nodes_table = res$nodes ,
									  ges_signature_vector = ges_signature_model_bad_match ,
									  pas_signature_vector = pas_signature_model_bad_match)
nodes_model_bad_match <- res$nodes
nodes_model_bad_match[ nodes_model_bad_match$type == "TF" , ]
filename <- file.path(reports.dir, paste0("blobs-",run_name,"-model-bad-match.pdf"))
plotBlobs( edges_table = res$edges , nodes_table = res$nodes , blob_filename = filename)

## TUMOR --
res <- generateNodesAndEdges(net=net_data,mrs_list=all_mrs)
res <- integrateExpressionAndActivity(edges_table = res$edges , nodes_table = res$nodes ,
									  ges_signature_vector = ges_signature_tumor_bad_match ,
									  pas_signature_vector = pas_signature_tumor_bad_match)
nodes_tumor_bad_match <- res$nodes
nodes_tumor_bad_match[ nodes_tumor_bad_match$type == "TF" , ]
filename <- file.path(reports.dir, paste0("blobs-",run_name,"-tumor-bad-match.pdf"))
plotBlobs( edges_table = res$edges , nodes_table = res$nodes , blob_filename = filename)






