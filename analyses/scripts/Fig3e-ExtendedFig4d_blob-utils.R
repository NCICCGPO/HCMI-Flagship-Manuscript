##
# Helper function to generate Blobs OncoLoop-style
# ------------------------------------------------

## Function to generate tables of nodes and edges from ARACNe network ----
generateNodesAndEdges <- function(net,mrs_list=NULL)
{
	if ( !is.null(mrs_list) )
	{
		# pruning networks of regulators by retaining only a subset of interest
		net <- net[ names(net) %in% mrs_list ]		
	}
	
	print_msg_info(">>> Interactome conversion to data frame (tibble)")
	edges <- transformInteractomeToTibble(net,regulon_type = "narnea")
	
	print_msg_info(">>> >>> Handling Edges ...")
	edges <- dplyr::rename(edges, source = TF)
	edges <- dplyr::rename(edges, target = target_gene_name)
	
	edges$id <- 1:nrow(edges)
	
	print_msg_info(">>> >>> Handling Nodes ...")
	nodes <- tibble( id = unique( c( edges$source , edges$target ) ) , 
					 type = "" , 
					 stringsAsFactors = FALSE )
	nodes$type <- ifelse( nodes$id %in% edges$target , "TARGET" , nodes$type )		
	nodes$type <- ifelse( nodes$id %in% edges$source , "TF" , nodes$type )
	table(nodes$type)
	
	return(list(nodes=nodes,edges=edges))
}

## Function to Addi GES and PAS to edges and nodes table ----
integrateExpressionAndActivity <- function(edges_table,nodes_table,ges_signature_vector,pas_signature_vector)
{
	print_msg_info(">>> Integration with Gene Expression Signatures and Protein Activity ...")
	library(scales)
	
	print_msg_info(">>> Adding gene expression scores ...")
	# ges_tbl <- tibble( gene = names(ges_signature_model) , expression = ges_signature_model )
	ges_tbl <- tibble( gene = names(ges_signature_vector) , expression = ges_signature_vector )
	
	# print_msg_info(">>> Scaling gene expression scores ...")
	# summary(ges_tbl$expression)
	# ges_tbl$expression <- rescale(ges_tbl$expression,to=c(-10,10))
	# summary(ges_tbl$expression)
	
	## Calc Expression function ----
	edges_and_nodes.tibble <- left_join( edges_table , ges_tbl , by = c("target"="gene") )
	
	# edges_and_nodes.tibble$tfmode <- rescale(edges_and_nodes.tibble$tfmode,to=c(-1,1))
	# edges_and_nodes.tibble <- edges_and_nodes.tibble %>% mutate(expression_final=tfmode*expression)
	
	print_msg_info(">>> Adding protein activity scores ...")		
	# pas_df <- tibble( protein = names(pas_signature_model) , activity = pas_signature_model )
	pas_df <- tibble( protein = names(pas_signature_vector) , activity = pas_signature_vector )
	
	edges_and_nodes.tibble <- left_join( edges_and_nodes.tibble , pas_df , by = c("source"="protein") )
	
	nodes_table <- left_join( nodes_table , edges_and_nodes.tibble %>% dplyr::select(source,activity) %>% distinct() ,
							  by = c("id"="source"))
	# nodes_table <- left_join( nodes_table , edges_and_nodes.tibble %>% dplyr::select(target,expression_final) %>% distinct() ,
	nodes_table <- left_join( nodes_table , edges_and_nodes.tibble %>% dplyr::select(target,expression) %>% distinct() ,
							  by = c("id"="target"))
	# nodes_table <- nodes_table %>% group_by(id,type,activity) %>% summarise(across(expression_final, list(mean)))
	nodes_table <- nodes_table %>% group_by(id,type,activity) %>% summarise(across(expression, list(mean)))
	# nodes_table <- nodes_table %>% dplyr::rename(expression=expression_final_1)
	nodes_table <- nodes_table %>% dplyr::rename(expression=expression_1)
	
	return(list(nodes=nodes_table,edges=edges_table))
	
}

## Visualize Networks ----
plotBlobs <- function(edges_table , nodes_table , blob_filename , N_for_rescale = 10 )
{
	print_msg_info(">>> Network visualization with <<network>>")

	## TO FIX (not to use select distinct but to make unique keys)
	nodes_table <- left_join( nodes_table , 
							  edges_table %>% dplyr::select(-source,-id) , 
							  by = c("id"="target") ) # %>% dplyr::distinct(id,type,activity,expression,tfmode,likelihood)
	nodes_table <- nodes_table %>% group_by(id,type,activity) %>% 
		summarise( across(c(expression,tfmode,likelihood), list(mean) , .names = "{.col}" ) )
	# nodes_table$expression_to_plot <- nodes_table$expression * nodes_table$tfmode * nodes_table$likelihood * 10
	nodes_table$expression_to_plot <- nodes_table$expression * nodes_table$tfmode * 10
	
	nodes_table <- nodes_table %>% arrange(type)
	
	nodes_table$expression_to_plot[ nodes_table$expression_to_plot > 10 ] <- 10
	nodes_table$expression_to_plot[ nodes_table$expression_to_plot < -10 ] <- -10	
	nodes_table$activity[ nodes_table$activity > 10 ] <- 10
	nodes_table$activity[ nodes_table$activity < -10 ] <- -10
	
	N_for_rescale.act <- N_for_rescale
	N_for_rescale.expr <- N_for_rescale
	
	require(igraph)
	require(network)
	require(scales)
	require(viridis)
	
	# edges_table$weights <- sqrt(sqrt(edges_table$likelihood)*edges_table$tfmode^2) * 100
	edges_table$weights <- 1
	g <- graph_from_data_frame(d = edges_table,vertices = nodes_table , directed = FALSE)
	
	V(g)$vertex.label <- names(V(g))
	V(g)$vertex.label <- ifelse( V(g)$vertex.label %in% nodes_table$id[ nodes_table$type == "TF"] , V(g)$vertex.label , "" )
	
	set.seed(666)
	# Basic chart
	# plot(g)
	# coords <- layout_nicely(g,dim = 2)
	# coords <- layout.fruchterman.reingold(g)
	# coords <- layout.lgl( g , area = 1e8 , cellsize = 1e4 )
	# coords <- layout.reingold.tilford(g, circular=F)
	# coords <- layout.fruchterman.reingold(g, params = list(niter=1000) )
	
	coords <- layout.fruchterman.reingold(graph = g,weights = E(g)$weights , params = list(niter=1000))
	# coords <- layout_nicely(graph = g,weights = E(g)$weights )
	# coords <- layout_(g, as_bipartite())
	# coords <- layout.lgl( g , area = 1e8 , cellsize = 1e4 )
	
	# minC <- rep(-Inf, vcount(g))
	# maxC <- rep(Inf, vcount(g))
	# minC[1] <- maxC[1] <- 0
	# coords <- layout_with_fr(g, minx=minC, maxx=maxC,
	# 	miny=minC, maxy=maxC)
	
	# color palette
	library(RColorBrewer)
	# E(g)$color <- "black"
	
	# V(g)$color <- ifelse(V(g)$type == "TF", "darkred", "darkorange")
	V(g)$size <- ifelse(V(g)$type == "TF", 15 , 3.5 )
	
	# V(g)$color <- ifelse(V(g)$type == "TF", "darkred", "darkorange")
	
	# # V(g)$score_ges <- ifelse( V(g)$score_ges > 12 , 12 , V(g)$score_ges )
	# # V(g)$score_ges <- ifelse( V(g)$score_ges < -12 , -12 , V(g)$score_ges )
	# V(g)$score_ges_rescaled <- rescale( V(g)$score_ges )
	# V(g)$score_pas <- ifelse( V(g)$score_pas > 9 , 9 , V(g)$score_pas )
	# V(g)$score_pas <- ifelse( V(g)$score_pas < -9 , -9 , V(g)$score_pas )
	# V(g)$score_pas_rescaled <- rescale( V(g)$score_pas )		
	
	rangeColorFunc <- function(x,n=N_for_rescale) { 
		res <- 1+n*(x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))
		res <- ifelse( res == 0 , 1 , res )
		res <- ifelse( res >= N_for_rescale , N_for_rescale , res )
		return(res)
	}
	
	E(g)$color <- "gray"
	# colr <- viridis(N_for_rescale.expr)
	# colr <- magma(N_for_rescale.expr)
	colr <- inferno(N_for_rescale.expr)
	
	# colr <- colorRampPalette( c("dodgerblue3","gray50","firebrick2") )(N_for_rescale.expr)
	# V(g)$color[ V(g)$type == "TARGET" ] <- colr[ round( rangeColorFunc( V(g)$expression[V(g)$type == "TARGET"] , N_for_rescale.expr ) ) ]
	# colr <- colorRampPalette( c("dodgerblue3","gray50","firebrick2") )(N_for_rescale.act)
	# V(g)$color[ V(g)$type == "TF" ] <- colr[ round( rangeColorFunc( V(g)$activity[V(g)$type == "TF"] , N_for_rescale.act ) ) ]
	
	# Apply the custom color function to the "TARGET" nodes based on "expression"
	V(g)$color[V(g)$type == "TARGET"] <- BlobsScoresBinning(V(g)$expression_to_plot[V(g)$type == "TARGET"] , color_palette = "PuOr" )$colors
	# Apply the custom color function to the "TF" nodes based on "activity"
	V(g)$color[V(g)$type == "TF"] <- BlobsScoresBinning(V(g)$activity[V(g)$type == "TF"])$colors
	
	# E(g)$weights <- 0
	# E(g)$tfmode <- 0
	# E(g)$likelihood <- 1
	V(g)$activity <- 1
	V(g)$expression <- 1
	
	
	# g <- g %>% set_edge_attr("color", value = "red")
	
	# plot
	# par(bg="grey13", mar=c(0,0,0,0))
	# par( bg = "white", mar=c(0,0,0,0))
	pdf( blob_filename )
	my_plot <- plot(g,
					rescale = TRUE,
					layout = coords ,
					# vertex.size = c(5,2),
					vertex.label=V(g)$vertex.label,
					vertex.label.cex = 0.5,
					vertex.label.family = "sans" ,
					vertex.label.color = "white",
					# vertex.frame.color = "transparent" ,
					vertex.frame.color = "gray75" ,
					vertex.frame.size = 5 ,
					
					rescale=FALSE,
					
					# edge.color = rep(c("red","pink"),5),           # Edge color
					# edge.color = "gray",           # Edge color
					edge.width = 0.5,
					# edge.width = seq(1,10),                        # Edge width, defaults to 1
					edge.arrow.size = 0.35,                           # Arrow size, defaults to 1
					edge.arrow.width = 0.15,                          # Arrow width, defaults to 1
					edge.lty = c("solid")                           # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
					#edge.curved=c(rep(0,5), rep(1,5))            # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
	)
	dev.off()
	
	# pdf( file.path(reports.dir,"blobs.pdf") )
	# 	print(my_plot)
	# dev.off()
	
	# par(mfrow=c(2,2), mar=c(1,1,1,1))
	# plot(g, layout=layout.sphere, main="sphere")
	# plot(g, layout=layout.circle, main="circle")
	# plot(g, layout=layout.random, main="random")
	# plot(g, layout=layout.fruchterman.reingold, main="fruchterman.reingold")		
	
	# library(ggplot2)
	# qplot(x=-10:10, y = 1, fill=-10:10, geom="tile") +
	# 	# scale_fill_manual(values = rgb(ddf$r, ddf$g, ddf$b)) +
	# 	# scale_fill_gradient2(low="dodgerblue3",mid="gray50",high="firebrick2") +
	# 	scale_fill_gradient2(low=colr[1],mid=colr[5],high=colr[10]) +
	# 	theme_void()+
	# 	theme(legend.position="none") 
	
}

BlobsScoresBinning <- function( vector , 
								   color_palette = "RdBu" ,
								   breaks = NULL ,
								   tags = NULL ) {
	
	# my_max <- round(max(mat),0)
	# set up cut-off values 
	# if ( my_max <= 30 && (is.null(breaks) || is.null(tags)) )
	# {
		breaks <- c(-10,-5,-3,-2,0,+1,+2,+3,+5,+10)		
		# specify interval/bin labels
		tags <- c("[-10-5)","[-5-3)", "[-3-2)", "[-2-1)", "[0-0)", "[+1+2)", "[+2+3)" , "[+3+5)" ,"[+5+10)")
		tags <- factor(tags)    
	# } else if (is.null(breaks) || is.null(tags)) {
	# 	breaks <- c(0,5,10,15,20,25,30,my_max+1)		
	# 	# specify interval/bin labels
	# 	tags <- c("[0-5)","[5-10)", "[10-15)", "[15-20)", "[20-25)", "[25-30)",paste0("[30-",my_max+1,")") )
	# 	tags <- factor(tags)    
	# }
	
	# x <- RColorBrewer::brewer.pal( length(breaks)-1 , "Set1" )
	color_map <- rev(RColorBrewer::brewer.pal( length(breaks)-1 , color_palette ))
	# x[1] <- "#FFFFFF"
	names(color_map) <- tags
	# f1 = colorRamp2( as.integer(tags) , x )
	
	# group_tags <- apply(mat, 2 , function(x) { 
		res <- cut(vector, 
				   breaks=breaks, 
				   include.lowest=TRUE, 
				   right=FALSE, 
				   labels=tags ,
				   # labels=FALSE ,
				   ordered_result=TRUE) ;
		# res <- factor(res,levels = tags,ordered = TRUE)
	# })  
	# rownames(group_tags) <- rownames(mat)
	names(res) <- names(vector)
	
	return(list(tags=res,colors=color_map[res],scores=vector))
}

# vector <- rep(c(-10,-7,3,7,20),3)
# names(vector) <- vector
# BlobsScoresBinning(vector)
