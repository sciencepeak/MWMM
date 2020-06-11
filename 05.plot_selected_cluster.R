rm(list = ls())
ptm <- proc.time()
options(stringsAsFactors = F)

library("reshape2")
library("magrittr")
library("igraph")
source(file.path(getwd(), "function_scripts_directory", "draw_cluster.R"))


# To draw a cluster network, these parameters need to be designated.
cancer_type <- "BRCA"
algorithm_name <- "hungarian_algorithm"
desired_cluster_ordinal <- 189


# ---------read the raw edge list of integrated mean value weight.------------
raw_input_file_path <- file.path(getwd(),
                                 "output_data_directory",
                                 "plain_text_files",
                                 paste("six_kinds_of_weights_files", cancer_type, sep = "."),
                                 "integrated_mean_value_weight.csv")

raw_edge_list_df <- read.csv(raw_input_file_path, check.names = F)
raw_edge_list_df <- raw_edge_list_df[raw_edge_list_df$mRNA != "SLC35E2", ]

# make the raw graph and corresponding adjacency matrix from the input edge list.
raw_graph <- graph.data.frame(raw_edge_list_df, directed = F)
raw_adjacency_matrix <- as_adjacency_matrix(raw_graph, type="both", names=TRUE, sparse=FALSE, attr="weight")


# ----------------read the communities of a algorithm in a cancer.-------------
path_to_input_file <- file.path(getwd(),
                                 "output_data_directory",
                                 "binary_object_files",
                                 paste("communities_files", cancer_type, sep = "."),
                                 paste("the_communities", algorithm_name, "rds", sep = "."))
a_communities <- readRDS(path_to_input_file)


# ------------extract the edge list denoting the graph that you want to draw.----
cluster_nodes <- communities(a_communities)[[desired_cluster_ordinal]]

cluster_miRNAs <- cluster_nodes[grepl("hsa-", cluster_nodes)]
cluster_mRNAs <- cluster_nodes[!grepl("hsa-", cluster_nodes)]

condition_miRNA <- raw_edge_list_df$microRNA %in% cluster_miRNAs
condition_mRNA <- raw_edge_list_df$mRNA %in% cluster_mRNAs
condition_both <- condition_miRNA & condition_mRNA

cluster_edge_list_df <- raw_edge_list_df[which(condition_both), ]


# ------------draw the graph to the hard disk.------------------------------
file_type <- "pdf"
output_graph_file_name <- paste("cluster_bipartite_graph", algorithm_name, desired_cluster_ordinal, file_type, sep = ".")
path_to_graph_output <- file.path(getwd(), "output_data_directory", "graph_files", output_graph_file_name)
draw_cluster(cluster_edge_list_df, raw_edge_list_df, path_to_graph_output)


# ----------------Also draw the traditional hierarchical graph---------------
# sort the edgelist dataframe descreasingly by weight:  (highest -> lowest)
sorted_edgelist_df <- raw_edge_list_df[(order(raw_edge_list_df$weight, decreasing = T)), ]
hierarchical_edge_list_df <- sorted_edgelist_df[1:38, ]

output_graph_file_name2 <- paste("cluster_bipartite_graph", "hierarchical_cluster", "38", file_type, sep = ".")
path_to_graph_output2 <- file.path(getwd(), "output_data_directory", "graph_files", output_graph_file_name2)
draw_cluster(hierarchical_edge_list_df, raw_edge_list_df, path_to_graph_output2)


proc.time() - ptm



# -------------garbage codes---------------------------------------------------
# by the way, in the igraph package,
# graph_from_incidence_matrix creates a bipartite igraph graph from an incidence matrix. Try this latter. layout = bipartite...
# it is worth trying later...
# http://kateto.net/network-visualization


