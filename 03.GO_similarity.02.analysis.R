rm(list=ls())
ptm <- proc.time()
options(stringsAsFactors = F)
library("magrittr")
library("tidyr")


# look up the input directory to get the cancer types.
TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])


for (cancer_type in cancer_types) {
    
    
    path_to_GO_similarity_text_file_directory <- file.path(getwd(),
                                                           "output_data_directory",
                                                           "plain_text_files",
                                                           paste("GO_similarity_distance", cancer_type, sep = "."))
    
    # if the output file is there, remove it,
    # so each run of this script doesn't append to the previous input.
    if (file.exists(file.path(path_to_GO_similarity_text_file_directory, "GO_similarity.summary.csv"))) {
        file.remove(file.path(path_to_GO_similarity_text_file_directory, "GO_similarity.summary.csv"))
    }
    
    # get algorithm names, ie, communities_rounds.
    communities_rounds <- dir(path_to_GO_similarity_text_file_directory) %>%
        strsplit(., "\\.") %>%
        sapply(., function(x) x[2][[1]])
    
    
    # read the calcuated GO similarity text files to a list for each cancer type.
    input_file_names <- paste("GO_similarity", communities_rounds, "csv", sep = ".")
    
    path_to_input_files_vector <- file.path(path_to_GO_similarity_text_file_directory,
                                            input_file_names)
    
    a_list_of_GO_similarity_df <- lapply(path_to_input_files_vector, read.csv, check.names = F)
    
    names(a_list_of_GO_similarity_df) <- communities_rounds
    
    
    accum_list <- vector(mode = "list", length = length(communities_rounds))
    for(i in names(a_list_of_GO_similarity_df)){
        
        current_GO_df <- a_list_of_GO_similarity_df[[i]]
        
        # Add a column that is the intra-cluster score minus inter cluster score.
        current_GO_df$score_difference <- current_GO_df$`in-cluster_similarity_score` - current_GO_df$`inter-cluster_similarity_score`
        
        # Calculate the the mean of each column grouped by the GO terms: BP, CC, MF
        wide_df <- aggregate(current_GO_df[, -c(1, 2)], list(current_GO_df$ontology), mean)
        colnames(wide_df) = c("GO_term", "intra-cluster", "inter-cluster", "difference")
        
        # organize the data frame from the wide format table to the long format table.
        long_df <- gather(data = wide_df,
                          key = "score_category",
                          value = score_value, "intra-cluster", "inter-cluster", difference)
        
        sorted_long_df <- long_df[order(long_df$GO_term), ]
        rownames(sorted_long_df) <- paste(sorted_long_df$GO_term,
                                          sorted_long_df$score_category, sep = "_")
        
        # only keep the average score values.
        concise_df <- sorted_long_df[, 3, drop = F]
        transposed_matrix <- t(concise_df)
        transposed_vector <- transposed_matrix[1, ]
        
        accum_list[[i]] <- transposed_vector
        cat(cancer_type, i, "has be summarized", sep = " ", "\n")
        print("-------------------------------")
    }
    
    all_algorithm_similarity_score_df <- do.call(rbind, accum_list)
    
    path_to_output_file <- file.path(path_to_GO_similarity_text_file_directory, "GO_similarity.summary.csv")
    write.csv(all_algorithm_similarity_score_df, file = path_to_output_file, quote = F, row.names = T)
    cat(cancer_type, "has done for GO similarity summary.", sep = " ", "\n")
    print("***********************************************************")

}

proc.time() - ptm

