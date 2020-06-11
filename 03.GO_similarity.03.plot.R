rm(list=ls())
options(stringsAsFactors = F)
library("ggplot2")
library("tidyr")


# look up the input directory to get the cancer types.
TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])

graph_output_path <- file.path(getwd(), "output_data_directory", "graph_files", "GO_similarity_scores_graphs")
dir.create(graph_output_path, showWarnings = F, recursive = T)

for (cancer_type in cancer_types) {
    
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
    

    path_to_GO_similarity_summary_file <- file.path(getwd(),
                                                    "output_data_directory",
                                                    "plain_text_files",
                                                    paste("GO_similarity_distance", cancer_type, sep = "."),
                                                    "GO_similarity.summary.csv")
    all_algorithm_similarity_score_df <- read.csv(path_to_GO_similarity_summary_file, check.names = F)
    colnames(all_algorithm_similarity_score_df)[1] <- "Clustering_algorithm"
    
    BP_term_wide_format_df <- all_algorithm_similarity_score_df[, c(1, 2, 3, 4)]
    CC_term_wide_format_df <- all_algorithm_similarity_score_df[, c(1, 5, 6, 7)]
    MF_term_wide_format_df <- all_algorithm_similarity_score_df[, c(1, 8, 9, 10)]
    
    biological_term_df_list <- list(BP_term_wide_format_df,
                                    CC_term_wide_format_df,
                                    MF_term_wide_format_df)
    names(biological_term_df_list) <- c("BP", "CC", "MF")
    
    for(term_name in names(biological_term_df_list)) {

        one_term_wide_format_df <- biological_term_df_list[[term_name]]
        
        one_term_long_format_df <- gather(data = one_term_wide_format_df,
                                          key = category,
                                          value = Average_GO_similarity_distance_score,
                                          colnames(one_term_wide_format_df)[2:4])
        one_term_long_format_df$category <- factor(one_term_long_format_df$category, levels = colnames(one_term_wide_format_df)[2:4])
        
        one_term_long_format_df$Clustering_algorithm <- factor(one_term_long_format_df$Clustering_algorithm, levels = communities_rounds)
        
        
        # Grouped Bar Plot
        ggplot(one_term_long_format_df,
               aes(fill=category, y=Average_GO_similarity_distance_score, x=Clustering_algorithm)) + 
            geom_bar(position="dodge", stat="identity") + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        file_types <- c("pdf", "wmf")
        
        for(file_type in file_types) {
            output_graph_file_name <- paste("GO_term_similarity_distance_score", cancer_type, term_name, file_type, sep = ".")
            
            ggsave(output_graph_file_name, plot = last_plot(),
                   width = 170, height = 170, units = "mm", device = file_type,
                   path = graph_output_path)
        }
    }
}

