rm(list=ls())
ptm <- proc.time()
options(stringsAsFactors = F)
library("magrittr")
library("igraph")


# We only run BRCA to compare six weights formula on one algorithm results.
cancer_type <- "BRCA"
algorithm_name <- "blossom_01"
# algorithm_name <- "hungarian_algorithm"

# Setup input
six_weight_formula_vector <- c("all_negative_value_weight",
                               "all_positive_value_weight",
                               "arithmetic_mean_value_weight",
                               "geometric_mean_value_weight",
                               "integrated_mean_value_weight",
                               "maximum_absolute_value_weight")
communities_output_dir <- file.path(getwd(),
                                    "output_data_directory",
                                    "binary_object_files",
                                    "six_weights",
                                    six_weight_formula_vector,
                                    paste("communities_files", cancer_type, sep = "."),
                                    paste("the_communities", algorithm_name, "rds", sep = "."))


# set up output for adjusted rand index (ARI) results.
ARI_comparison_out_put_dir <- file.path("output_data_directory", "plain_text_files",
                                        "six_weights", "ajusted_rand_index_results")
dir.create(ARI_comparison_out_put_dir, showWarnings = F, recursive = T)
ARI_output_file_name <- paste("six_weight_ARI", algorithm_name, cancer_type, "csv", sep = ".")
output_dir_and_name <- file.path(ARI_comparison_out_put_dir, ARI_output_file_name)


# Read the input rds files that contain communities structure in igraph.
# Read the rds into a named list, in which name is the six formula names.
communities_list <- lapply(communities_output_dir, readRDS)
names(communities_list) <- six_weight_formula_vector


# calculate adjusted rand index (ARI) of every two communities structures.
combination_matrix <- combn(six_weight_formula_vector, 2)


calculate_adjusted_rand_index <- function(a_column_of_a_matrix) {
    
    weight_formula_A <- a_column_of_a_matrix[1]
    weight_formula_B <- a_column_of_a_matrix[2]
    
    communities_A <- communities_list[[weight_formula_A]]
    communities_B <- communities_list[[weight_formula_B]]
    # call the compare function in igraph package to compare the distance between communities.
    # the method is adjusted rand index (ARI)
    ARI <- compare(communities_A, communities_B, method = "adjusted.rand")
    comparison_vector <- c(weight_formula_A, weight_formula_B, ARI)
    return(comparison_vector)
}

# got the edge list data frame that have the every two communities structure ARI as weights.
edgelist_of_weight_formulas_df <- apply(combination_matrix, 2, calculate_adjusted_rand_index) %>%
    t %>%
    as.data.frame %>%
    set_colnames(., c("from", "to", "weight"))
edgelist_of_weight_formulas_df[, "weight"] <- round(as.numeric(edgelist_of_weight_formulas_df$weight), digits = 3)


# Write results to hard disk.
edgelist_of_weight_formulas_df3 <- edgelist_of_weight_formulas_df
colnames(edgelist_of_weight_formulas_df3) <- c("from", "to", "ARI")
write.csv(edgelist_of_weight_formulas_df3, output_dir_and_name, row.names = F, quote = F)


# following is the block to draw the dendrogram to see the dissimilarity of formulas' consequences.
edgelist_of_weight_formulas_df2 <- edgelist_of_weight_formulas_df
edgelist_of_weight_formulas_df2[, "weight"] <- scale (1 / edgelist_of_weight_formulas_df$weight)

# make the raw graph and corresponding adjacency matrix from the input edge list.
raw_graph <- graph.data.frame(edgelist_of_weight_formulas_df2, directed = F)
raw_adjacency_matrix <- as_adjacency_matrix(raw_graph, type="both", names=TRUE, sparse=FALSE, attr="weight")

# Convert the adjacency matrix to the dist object.
# https://stackoverflow.com/questions/17875733/how-to-convert-a-symmetric-matrix-into-dist-object
dist_obj <- as.dist(raw_adjacency_matrix)

hclust_obj <- hclust(dist_obj, method = "average")
plot(hclust_obj)



# ----------garbage code------------------------
# g <- make_graph("Zachary")
# sg <- cluster_spinglass(g)
# le <- cluster_leading_eigen(g)
# compare(sg, le, method="adjusted.rand")
# compare(membership(sg), membership(le))
# ## Zachary's karate club
# g <- make_graph("Zachary")
# 
# imc <- cluster_infomap(g)
# membership(imc)
# communities(imc)
# 
# 
# library("clues")
# set.seed(1)
# x <- sample(x = rep(1:3, 4), 12)
# set.seed(2)
# y <- sample(x = rep(1:3, 4), 12)
# 
# 
# 
# adjustedRand(x, y)
# library("mclust")
# adjustedRandIndex(x, y)
# 
# compare(x, y, method = "adjusted.rand")
# compare(communities_list[["hungarian_algorithm"]], communities_list[["hungarian_algorithm"]], method = "adjusted.rand")

# # make the communities from a graph
# # if you have a graph, the igraph function clusters can make a communities for you.
# list_format_cluster <- clusters(raw_graph)
# 
# # clusters() can only make a list-format cluster, make_clusters() makes the real communities object.
# communities_format_cluster <- make_clusters(graph = raw_graph,
#                                             membership = list_format_cluster$membership)
# 
# membership(communities_format_cluster)
