rm(list=ls())
options(stringsAsFactors = F)
ptm <- proc.time()
args<-commandArgs(TRUE)
library("magrittr")

# -----------------------------------------------------------------------------
source(file.path(getwd(), "function_scripts_directory", "calculate_enrichment_from_a_list_of_cluster.R"))
source(file.path(getwd(), "function_scripts_directory", "GO_similarity_functions.R"))
#------------------source functions below---------------------------------------
conversion_table_path.TCGA <- file.path(getwd(), "function_scripts_directory", "mRNA_names_and_Entrez_IDs_dictionary.TCGA.csv")

TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])

for (cancer_type in cancer_types) {
    
    output_directory <- file.path(getwd(),
                                  "output_data_directory",
                                  "binary_object_files",
                                  paste("enrichment_analysis", cancer_type, sep = "."))
    dir.create(output_directory, showWarnings = F, recursive = T)
    
    
    # get algorithm names, ie, communities_rounds.
    path_to_communities_input_directory <- file.path(getwd(),
                                                     "output_data_directory",
                                                     "binary_object_files",
                                                     paste("communities_files", cancer_type, sep = "."))
    
    
    # get the names of different algorithms.
    communities_rounds <- dir(path_to_communities_input_directory) %>%
        strsplit(., "\\.") %>%
        sapply(., function(x) x[2][[1]])
    
    
    input_file_names <- paste("the_communities", communities_rounds, "rds", sep = ".")
    path_to_input_files <- file.path(getwd(),
                                     "output_data_directory",
                                     "binary_object_files",
                                     paste("communities_files", cancer_type, sep = "."),
                                     input_file_names)
    
    a_list_of_the_communities <- lapply(path_to_input_files, readRDS)
    names(a_list_of_the_communities) <- communities_rounds
    
    
    # -------------convert clusters in communities object to a list of list of clusters
    a_list_of_cluster_list <- vector(mode = "list", length = length(communities_rounds))
    
    for (i in seq_along(communities_rounds)) {
        path_to_input_file <- path_to_input_files[i]
        current_communities <- readRDS(path_to_input_file)
        
        current_cluster_list <- convert_a_communities_to_a_list_of_mRNA_only_clusters(current_communities, conversion_table_path.TCGA)
        
        # filter out the clusters that only have miRNAs.
        current_cluster_list <- current_cluster_list[!sapply(current_cluster_list, is.list)]
        
        a_list_of_cluster_list[[i]] <- current_cluster_list
        print(i)
    }
    names(a_list_of_cluster_list) <- communities_rounds

    
    for (algorithm_name in names(a_list_of_cluster_list)) {
        this_cluster_list <- a_list_of_cluster_list[[algorithm_name]]
        
        current_enriched_terms_list <- calculate_enrichment_from_a_cluster_list(this_cluster_list)
        output_name <- paste("biological_term_list", algorithm_name, "rds", sep = ".")
        output_object_path <- file.path(output_directory, output_name)
        saveRDS(current_enriched_terms_list, file = output_object_path)
        
        cat(cancer_type, algorithm_name, "is done", sep = " ", "\n")
    }
}

proc.time() - ptm




# ----------------------garbage codes---------------------------------
# all_indicator_list <- foreach(algorithm_name = names(a_list_of_cluster_list)) %dopar% {
# the_cluster_list <- a_list_of_cluster_list$blossom_02

# path_to_output_file <- file.path(getwd(), "output_data_directory", "plain_text_files", "GO_similarity_distance.BRCA", output_file_names[i])
# write.csv(GO_similarity_df, file = path_to_output_file, row.names = F, quote = F)

# cancer_types <- rev(cancer_types)
# conversion_table_path.miRMAP_bicluster <- file.path(getwd(), "function_scripts_directory", "mRNA_names_and_Entrez_IDs_dictionary.miRMAP_bicluster.csv")


# if (communities_rounds[i] == "miRMAP_bicluster") {
#     current_cluster_list <- convert_a_communities_to_a_list_of_mRNA_only_clusters(current_communities, conversion_table_path.miRMAP_bicluster)
# } else {
#     current_cluster_list <- convert_a_communities_to_a_list_of_mRNA_only_clusters(current_communities, conversion_table_path.TCGA)
# }
