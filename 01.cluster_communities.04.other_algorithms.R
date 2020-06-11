rm(list=ls())
ptm <- proc.time()
options(stringsAsFactors = F)
library("igraph")
library("magrittr")
library("foreach")
library("doParallel") # the parallel package will be loaded by doParallel


# look up the input directory to get the cancer types.
TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])

# -----------------------------------------------------------------------------
# --------------------- Set up parallel running below: ------------------------


# parameters to change on different system and machines:
# extra_libraries_path_directory <- "/net/home/lding/R/x86_64-redhat-linux-gnu-library/3.4/"

# my packages are also installed on the following directories.
# .libPaths(c(.libPaths(), extra_libraries_path_directory))

#setup parallel backend to use many processors, try all the available cores?
available_cores_number <- detectCores()
if (available_cores_number > length(cancer_types) ) {
    available_cores_number <- length(cancer_types)
}

parallel_cluster <- makeCluster(available_cores_number)
registerDoParallel(parallel_cluster)

# load packages on each worker, so that all the threads can use packages by parLapply()
# http://stackoverflow.com/questions/17196261/understanding-the-differences-between-mclapply-and-parlapply-in-r
clusterEvalQ(parallel_cluster, library("igraph"))


final_result <- foreach(cancer_type = cancer_types) %dopar% {
    
# for (cancer_type in cancer_types) {
    
    print(cancer_type)
    
    raw_input_file_path <- file.path(getwd(),
                                     "output_data_directory",
                                     "plain_text_files",
                                     paste("six_kinds_of_weights_files", cancer_type, sep = "."),
                                     "integrated_mean_value_weight.csv")
    
    raw_edge_list_df <- read.csv(raw_input_file_path, check.names = F)[, 1:3]
    raw_edge_list_df <- raw_edge_list_df[raw_edge_list_df$mRNA != "SLC35E2", ]
    
    # make the raw graph and corresponding adjacency matrix from the input edge list.
    raw_graph <- graph.data.frame(raw_edge_list_df, directed = F)
    raw_adjacency_matrix <- as_adjacency_matrix(raw_graph, type="both", names=TRUE, sparse=FALSE, attr="weight")
    
    # components(raw_graph)
    
    
    # create output directory
    output_directory_path <- file.path(getwd(),
                                  "output_data_directory",
                                  "binary_object_files",
                                  paste("communities_files", cancer_type, sep = "."))
    dir.create(output_directory_path, showWarnings = F, recursive = T)
    
    
    # Functions to deal with the result of network community detection
    # http://igraph.org/r/doc/communities.html
    # 
    # What are the differences between community detection algorithms in igraph?
    # https://stackoverflow.com/questions/9471906/what-are-the-differences-between-community-detection-algorithms-in-igraph
    
    
    
    # -------------------------Louvain cluster algorithm------------------------
    # The louvain clustering algorithm returns an object called communities in R.
    # REF [Functions to deal with the result of network community detection](http://igraph.org/r/doc/communities.html)
    the_communities.louvain_algorithm <- cluster_louvain(raw_graph)
    
    path_to_output_file <- file.path(output_directory_path, "the_communities.louvain_algorithm.rds")
    
    saveRDS(the_communities.louvain_algorithm, file = path_to_output_file)
    
    
    
    # -------------------------fast_greedy------------------------------------
    the_communities.fast_greedy.trial <- try(cluster_fast_greedy(raw_graph, weights = E(raw_graph)$weight))
    
    
    if (is.null(attr(the_communities.fast_greedy.trial, "condition")$message)) {
        the_communities.fast_greedy <- the_communities.fast_greedy.trial
        print("fast greedy is done without any error")
    } else {
        simple_graph <- simplify(raw_graph)
        the_communities.fast_greedy <- cluster_fast_greedy(simple_graph, weights = E(simple_graph)$weight)
        print( "fast_greedy is done by simplifying the graph.")
    }
    
    path_to_output_file <- file.path(output_directory_path, "the_communities.fast_greedy.rds")
    
    saveRDS(the_communities.fast_greedy, file = path_to_output_file)
    
    
    
    # --------------------------walktrap----------------------------------------
    the_communities.walktrap_algorithm <- cluster_walktrap(raw_graph, weights = E(raw_graph)$weight)
    
    path_to_output_file <- file.path(output_directory_path, "the_communities.walktrap_algorithm.rds")
    
    saveRDS(the_communities.walktrap_algorithm, file = path_to_output_file)
    
    
    
    # ------------------------leading_eigen-------------------------------------
    the_communities.leading_eigen <- cluster_leading_eigen(raw_graph, weights = E(raw_graph)$weight, options=list(maxiter=100000))
    
    path_to_output_file <- file.path(output_directory_path, "the_communities.leading_eigen.rds")
    
    saveRDS(the_communities.leading_eigen, file = path_to_output_file)
    
    
    
    # -----------------------label-propagation----------------------------------
    the_communities.label_propagation <- cluster_label_prop(raw_graph, weights = E(raw_graph)$weight)
    
    path_to_output_file <- file.path(output_directory_path, "the_communities.label_propagation.rds")
    
    saveRDS(the_communities.label_propagation, file = path_to_output_file)
    
    
    
    # ------------------------cluster_edge_betweenness-------------------------
    the_communities.edge_betweenness <- cluster_edge_betweenness(raw_graph, weights = E(raw_graph)$weight, directed = FALSE)
    
    path_to_output_file <- file.path(output_directory_path, "the_communities.edge_betweenness.rds")
    
    saveRDS(the_communities.edge_betweenness, file = path_to_output_file)
    
    cat(cancer_type, "are all done.", sep = " ", "\n")

    status <- "successful"
    status
    
}
stopCluster(parallel_cluster)

print("Everything is done!")

proc.time() - ptm