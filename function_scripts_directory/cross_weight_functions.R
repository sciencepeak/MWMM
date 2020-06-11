library("igraph")
library("magrittr")
library("foreach")
library("doParallel") # the parallel package will be loaded by doParallel


calculate_cross_weight_from_a_communities <- function(a_communities, reference_edge_list_df) {
    
    # parameters to change on different system and machines:
    # extra_libraries_path_directory <- "/net/home/lding/R/x86_64-redhat-linux-gnu-library/3.4/"
    
    # my packages are also installed on the following directories.
    # .libPaths(c(.libPaths(), extra_libraries_path_directory))
    
    #setup parallel backend to use many processors, use all the available cores
    available_cores_number <- detectCores()
    parallel_cluster <- makeCluster(available_cores_number)
    registerDoParallel(parallel_cluster)
    
    # load packages on each worker, so that all the threads can use packages by parLapply()
    # http://stackoverflow.com/questions/17196261/understanding-the-differences-between-mclapply-and-parlapply-in-r
    clusterEvalQ(parallel_cluster, library("igraph"))
    clusterEvalQ(parallel_cluster, library("magrittr"))

	# ---------------------------------------------------------------
	# send the variable to the environment on the nodes.
	clusterExport(cl = parallel_cluster,
	              c(deparse(substitute(a_communities)),
	                deparse(substitute(reference_edge_list_df))))

	# --------------------------------------------------------------------
	
	# get unique miRNA names
	cluster_ordinal_vector = seq_along(communities(a_communities))
	
    
    # If there are less than two clusters, we can not calculate cross weight.
    if (length(cluster_ordinal_vector) < 2 ) {
        cross_weight_edge_list_df <- data.frame(from = numeric(0),
                                                to = numeric(0),
                                                weight = numeric(0))
        
        return(cross_weight_edge_list_df)
    }
    
    
    
	cluster_pair_matrix <- combn(cluster_ordinal_vector, 2)
	
	# This function take the column of combination of the cluster ordinals.
	# Return the cross weights of each cluster pair combination.
	calculate_cross_weight_of_cluster_pairs <- function(a_column_of_cluster_pair_matrix) {
	    
	    # a_column_of_cluster_pair_matrix <- cluster_pair_matrix[, 1]
	    cluster_01 <- a_column_of_cluster_pair_matrix[1]
	    cluster_02 <- a_column_of_cluster_pair_matrix[2]
	    
	    # a_cluster_nodes is a character vector of mRNA and mRNA nodes
	    # a_cluster_nodes also represent a community/cluster in the communities object.
	    cluster_01_nodes <- communities(a_communities)[[cluster_01]]
	    cluster_02_nodes <- communities(a_communities)[[cluster_02]]
	    
	    # Distinguish the miRNA and mRNA nodes by prefix of miRNA.
	    cluster_01_miRNAs <- cluster_01_nodes[grepl("hsa-", cluster_01_nodes)]
	    cluster_02_miRNAs <- cluster_02_nodes[grepl("hsa-", cluster_02_nodes)]
	    
	    cluster_01_mRNAs <- cluster_01_nodes[!grepl("hsa-", cluster_01_nodes)]
	    cluster_02_mRNAs <- cluster_02_nodes[!grepl("hsa-", cluster_02_nodes)]
	    
	    
	    # get the logic index to subset the involved miRNA and mRNAs
	    condition_A_miRNA <- reference_edge_list_df$microRNA %in% cluster_01_miRNAs
	    condition_A_mRNA <- reference_edge_list_df$mRNA %in% cluster_02_mRNAs
	    condition_A_both <- condition_A_miRNA & condition_A_mRNA
	    
	    # get the logic index to subset the involved miRNA and mRNAs
	    condition_B_miRNA <- reference_edge_list_df$microRNA %in% cluster_02_miRNAs
	    condition_B_mRNA <- reference_edge_list_df$mRNA %in% cluster_01_mRNAs
	    condition_B_both <- condition_B_miRNA & condition_B_mRNA
	    
	    # condition A is the first cluster miRNA connect to the likely mRNAs in the second cluster.
	    # condition B is the second cluster miRNA connect to the likely mRNAs in the first cluster.
	    # OR operation is to select the edges that meet condition A or B.
	    # ie., to include both condition A and B edges.
	    all_condition_both <- condition_A_both | condition_B_both
	    
	    if (sum(all_condition_both) == 0) {
	        cross_weight_of_two_clusters <- 0
	    } else {
	        cross_weight_of_two_clusters <- reference_edge_list_df %>%
	            .[which(all_condition_both), ] %>%
	            .[, "weight"] %>%
	            sum
	        cross_weight_of_two_clusters <- cross_weight_of_two_clusters / (length(cluster_01_nodes) + length(cluster_02_nodes))
	    }
	    
	    return(cross_weight_of_two_clusters)
	}
	
	# in parallel apply the function to each combination.
	paired_cross_weight_vector <- parApply(cl = parallel_cluster,
	                                       cluster_pair_matrix,
	                                       2,
	                                       calculate_cross_weight_of_cluster_pairs)

	stopCluster(parallel_cluster)
	
	cross_weight_matrix <- rbind(cluster_pair_matrix, paired_cross_weight_vector)
	
	rownames(cross_weight_matrix) <- c("from", "to", "weight")
	
	non_zero_cross_weight_matrix <- cross_weight_matrix[, paired_cross_weight_vector != 0]
	
	cross_weight_edge_list_df <- as.data.frame(t(non_zero_cross_weight_matrix))
	
	return(cross_weight_edge_list_df)
}
