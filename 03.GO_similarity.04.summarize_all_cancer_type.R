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

path_to_GO_similarity_text_file_directory <- file.path(getwd(),
                                                       "output_data_directory",
                                                       "plain_text_files",
                                                       paste("GO_similarity_distance", cancer_types, sep = "."))

path_to_summary_files <- file.path(path_to_GO_similarity_text_file_directory, "GO_similarity.summary.csv")

a_list_of_GO_similarity_df <- lapply(path_to_summary_files, read.csv, row.names = 1, check.names = F)
names(a_list_of_GO_similarity_df) <- cancer_types


# find the top difference algorithms.
accum_list <- list()
for (cancer_type in cancer_types) {
    
    current_df <- a_list_of_GO_similarity_df[[cancer_type]]
    
    max_BP_difference_algorithm <- rownames(current_df)[which.max(current_df$BP_difference)]
    max_CC_difference_algorithm <- rownames(current_df)[which.max(current_df$CC_difference)]
    max_MF_difference_algorithm <- rownames(current_df)[which.max(current_df$MF_difference)]

    current_summary_vector <- c(cancer_type, max_BP_difference_algorithm,
                                max_CC_difference_algorithm, max_MF_difference_algorithm)
    names(current_summary_vector) <- c("cancer_type", "top_BP_difference", "top_CC_difference", "top_MF_difference")
    accum_list[[cancer_type]] <- current_summary_vector
}

all_cancer_df <- as.data.frame(do.call(rbind, accum_list))
new_all_cancer_df <- all_cancer_df
for (i in seq_len(nrow(all_cancer_df))) {
    current_row <- all_cancer_df[i, ]
    current_row[grepl("blossom|hungarian_algorithm", current_row)] <- "MWMM"
    new_all_cancer_df[i, ] <- current_row
}
rownames(new_all_cancer_df) <- seq_len(nrow(new_all_cancer_df))
write.csv(new_all_cancer_df,
          file = file.path(getwd(),
                           "output_data_directory",
                           "plain_text_files",
                           "all_cancer_type_GO_summary.csv"),
          row.names = F)

proc.time() - ptm


# -------------------garbage code--------------------------------------------
# find the top n difference algorithms.
# cutoff_n <- 3
# accum_list <- list()
# for (cancer_type in cancer_types) {
#     
#     current_df <- a_list_of_GO_similarity_df[[cancer_type]]
#     
#     large_BP_difference_algorithm <- current_df$BP_difference %>%
#         set_names(., rownames(current_df)) %>%
#         sort(., decreasing = T) %>%
#         head(., n = cutoff_n) %>%
#         names(.) %>%
#         paste(., collapse = ";")
#     
#     large_CC_difference_algorithm <- current_df$CC_difference %>%
#         set_names(., rownames(current_df)) %>%
#         sort(., decreasing = T) %>%
#         head(., n = cutoff_n) %>%
#         names(.) %>%
#         paste(., collapse = ";")
#     
#     large_MF_difference_algorithm <- current_df$MF_difference %>%
#         set_names(., rownames(current_df)) %>%
#         sort(., decreasing = T) %>%
#         head(., n = cutoff_n) %>%
#         names(.) %>%
#         paste(., collapse = ";")
#     
# 
#     current_summary_vector <- c(cancer_type,
#                                 large_BP_difference_algorithm,
#                                 large_CC_difference_algorithm,
#                                 large_MF_difference_algorithm)
#     
#     names(current_summary_vector) <- c("cancer_type",
#                                        paste("top", cutoff_n, "BP", "difference", sep = "_"),
#                                        paste("top", cutoff_n, "CC", "difference", sep = "_"),
#                                        paste("top", cutoff_n, "MF", "difference", sep = "_"))
#     accum_list[[cancer_type]] <- current_summary_vector
# }
# 
# all_cancer_df <- as.data.frame(do.call(rbind, accum_list))
# new_all_cancer_df <- all_cancer_df
# for (i in seq_len(nrow(all_cancer_df))) {
#     current_row <- all_cancer_df[i, ]
#     new_all_cancer_df[i, ] <- gsub("blossom|hungarian_algorithm", "MWMM", current_row)
# }
# rownames(new_all_cancer_df) <- seq_len(nrow(new_all_cancer_df))
# write.csv(new_all_cancer_df,
#           file = file.path(getwd(),
#                            "output_data_directory",
#                            "plain_text_files",
#                            "all_cancer_type_GO_summary_top_n.csv"),
#           row.names = F)


# communities_rounds <- c("MWMM", "fast_greedy", "leading_eigen", "edge_betweenness", "label_propagation", "louvain_algorithm", "walktrap_algorithm")
# refer to https://stackoverflow.com/questions/18774632/how-to-produce-stacked-bars-within-grouped-barchart-in-r



