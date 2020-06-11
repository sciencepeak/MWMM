rm(list=ls())
ptm <- proc.time()
options(stringsAsFactors = F)
library("igraph")
library("magrittr")

#------------------source functions below---------------------------------------
source(file.path(getwd(), "function_scripts_directory", "inner_and_outer_weight_functions.R"))

# -----------------input file setup.----------------------------------------------
# look up the input directory to get the cancer types.
TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])


for (cancer_type in cancer_types) {
    # cancer_type <- "BRCA"

    # get algorithm names, ie, communities_rounds.
    path_to_communities_input_directory <- file.path(getwd(),
                                                     "output_data_directory",
                                                     "binary_object_files",
                                                     paste("communities_files", cancer_type, sep = "."))
    
    
    # get the names of different algorithms.
    communities_rounds <- dir(path_to_communities_input_directory) %>%
        strsplit(., "\\.") %>%
        sapply(., function(x) x[2][[1]])
    
    # -----------------------------------------------------------------------
    
    # read the edge list for the current cancer type.
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
    
    # read the communities file for each algorithm, namely communities_round.
    communities_input_file_names <- paste("the_communities", communities_rounds, "rds", sep = ".")
    path_to_communities_input_files <- file.path(path_to_communities_input_directory, communities_input_file_names)
    a_list_of_the_communities <- lapply(path_to_communities_input_files, readRDS)
    names(a_list_of_the_communities) <- communities_rounds
    
    # prepare the output files for inner and outer weight for each algorithm.
    path_to_communities_output_directory <- file.path(getwd(),
                                                      "output_data_directory",
                                                      "plain_text_files",
                                                      paste("inner_and_outer_weight", cancer_type, sep = "."))
    dir.create(path_to_communities_output_directory, showWarnings = F, recursive = T)
    communities_output_file_names <- paste("inner_and_outer_weight", communities_rounds, "csv", sep = ".")
    path_to_weights_output_files <- file.path(path_to_communities_output_directory, communities_output_file_names)
    
    for (i in seq_along(communities_rounds)) {
        current_communities <- readRDS(path_to_communities_input_files[i])
        inner_and_outer_weight_df <- get_inner_outer_weights_and_conditions(current_communities, raw_edge_list_df)
        write.csv(inner_and_outer_weight_df, file = path_to_weights_output_files[i], quote = F, row.names = F)
    }
    
    # break
    cat(cancer_type, "has done!", sep = " ", "\n")
    print("*****************************************************************")
}
