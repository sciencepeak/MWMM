rm(list=ls())
options(stringsAsFactors = F)

ptm <- proc.time()

library("magrittr")
library("clusterProfiler")
library("ggplot2")

# look up the input directory to get the cancer types.
TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])


output_graph_directory <- file.path(getwd(), "output_data_directory", "graph_files", "cluster_enrichment_graphs")
dir.create(output_graph_directory, showWarnings = F, recursive = T)

for (cancer_type in cancer_types) {
    cancer_type <- "BRCA"
    
    # get algorithm names, ie, communities_rounds.
    path_to_communities_input_directory <- file.path(getwd(),
                                                     "output_data_directory",
                                                     "binary_object_files",
                                                     paste("communities_files", cancer_type, sep = "."))
    
    
    # get the names of different algorithms.
    communities_rounds <- dir(path_to_communities_input_directory) %>%
        strsplit(., "\\.") %>%
        sapply(., function(x) x[2][[1]])
    
    our_communities_rounds <- c("hungarian_algorithm", communities_rounds[grepl("blossom", communities_rounds)])
    other_communities_rounds <- c("fast_greedy",
                                  "leading_eigen",
                                  "edge_betweenness",
                                  "label_propagation",
                                  "louvain_algorithm",
                                  "walktrap_algorithm")
    
    communities_rounds <- c(other_communities_rounds, our_communities_rounds)
    
    
    input_file_names <- paste("biological_term_list", communities_rounds, "rds", sep = ".")
    path_to_input_files <- file.path(getwd(),
                                     "output_data_directory",
                                     "binary_object_files",
                                     paste("enrichment_analysis", cancer_type, sep = "."),
                                     input_file_names)
    
    a_list_of_enrichment_list <- lapply(path_to_input_files, readRDS)
    names(a_list_of_enrichment_list) <- communities_rounds
    
    
    for (i in names(a_list_of_enrichment_list)) {
        i <- "blossom_01"
        
        an_enrichment_result_list <- a_list_of_enrichment_list[[i]]
        
        #  an enrichment result of a clustering algorithm contains
        # "ck", "do", "ncg", "kk", "ego_BP", "ego_CC", "ego_MF"
        for (j in names(an_enrichment_result_list)) {
            
            j <- "kk"
            a_term_object <- an_enrichment_result_list[[j]]
            if (is.null(a_term_object)) {
                next
            } else {
                output_file_name <- paste("enrichment", cancer_type, i, j, "pdf", sep = ".")
                path_to_graph_output <- file.path(output_graph_directory, output_file_name)
                
                # The dotplot is calling ggplot2. If you understand Chinese, please
                # refer to https://guangchuangyu.github.io/cn/2017/07/clusterprofiler-dotplot/
                
                if (j == "kk") {
                    y_axis_label <- "Enriched pathways"
                } else {
                    y_axis_label <- "Enriched terms"
                }
                
                p <- dotplot(a_term_object)
                q <- p + labs(y=y_axis_label, x = "Clusters")
                # ggsave(path_to_graph_output, q, width = 170, height = 170, units = "mm", scale = 2.5)
            }
        }
    }
    break
}

temp_df <- a_term_object@compareClusterResult

proc.time() - ptm

