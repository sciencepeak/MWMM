rm(list=ls())
ptm <- proc.time()
options(stringsAsFactors = F)
library("magrittr")
library("igraph")
library("formattable") # This package provides percent function.

# This script compares the inner and outer weights of results from different algorithms.


cancer_type <- "BRCA"

# ------------------------------------------------------------------------
# setup input directories and files.

# We need to read the communities files so that we know the sizes of each communities,
# namely, how many nodes in each cluster in each communities.
communities_dir <- file.path(getwd(), "output_data_directory", "binary_object_files",
                             paste("communities_files", cancer_type, sep = "."))
communities_paths <- dir(communities_dir, full.names = T)

# We also need to read the inner weight and outer weight files to get the average.
weights_dir <- file.path(getwd(), "output_data_directory", "plain_text_files",
                         paste("inner_and_outer_weight", cancer_type, sep = "."))
weights_paths <- file.path(weights_dir, dir(weights_dir))

# setup the output directory and files.
output_dir <- file.path(getwd(),
                        "output_data_directory",
                        "plain_text_files",
                        "six_weights",
                        "inner_outer_weights_across_algorithms",
                        cancer_type)
dir.create(output_dir, showWarnings = F, recursive = T)
output_file <- "inner_outer_weights_assessment_across_algorithms.csv"
output_path <- file.path(output_dir, output_file)
# ----------------------------------------------------------------------------



# ----------------------------------------------------------------------------
# retrieve the algorithm names
algorithm_names <- dir(communities_dir) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[2])

# Read the communities structures and previously calculated inner and outer weights.
communities_list <- lapply(communities_paths, readRDS)
names(communities_list) <- algorithm_names

weights_list <- lapply(weights_paths, read.csv)
names(weights_list) <- algorithm_names
# ------------------------------------------------------------------------



#-------------------------------------------------------------------------
# calculation starts from here.
accum_list <- list()
for (current_algorithm_name in algorithm_names) {
    
    print(current_algorithm_name)

    current_weights_df <- weights_list[[current_algorithm_name]]
    current_communities <- communities_list[[current_algorithm_name]]
    
    # larger clusters have larger inner weights,
    # so they needs to be averages by node number of each cluster.
    cluster_sizes <- sizes(current_communities)
    normalized_IW <- current_weights_df$IW / cluster_sizes
    normalized_E2MROW <- current_weights_df$E2MROW / cluster_sizes
    normalized_R2MEOW <- current_weights_df$R2MEOW / cluster_sizes
    
    # to present the inner and outer weight, we calculate mean value.
    current_IW <- mean(normalized_IW)
    current_E2MROW <- mean(normalized_E2MROW)
    current_R2MEOW <- mean(normalized_R2MEOW)
    
    # Condition one specifies IW > E2MROW.
    # Condition two specifies 2 × IW > E2MROW + R2MEOW
    # condition_01 <- weighted_avg_IW > weighted_avg_E2MROW
    # condition_02 <- 2 * weighted_avg_IW > weighted_avg_E2MROW + weighted_avg_R2MEOW
    condition_01 <- normalized_IW > normalized_E2MROW
    condition_02 <- 2 * normalized_IW > normalized_E2MROW + normalized_R2MEOW

    # see how much percent of clusters that satisfy condition 1 and 2.
    condition_01_true_percent <- sum(condition_01)/length(condition_01)
    condition_02_true_percent <- sum(condition_02)/length(condition_02)
    

    current_vector <- c(current_IW, current_E2MROW, current_R2MEOW,
                        condition_01_true_percent, condition_02_true_percent)
    names(current_vector) <- c("average_IW", "average_E2MROW", "average_R2MEOW",
                               "condition_01_true_percent", "condition_02_true_percent")
    
    accum_list[[current_algorithm_name]] <- current_vector
    print("---------------------")
}

inner_outer_weights_summary_df <- do.call(rbind, accum_list)

# Format the numeric values in the data frame.
inner_outer_weights_summary_df[, 1] <- round(inner_outer_weights_summary_df[, 1], digits = 3)
inner_outer_weights_summary_df[, 2] <- round(inner_outer_weights_summary_df[, 2], digits = 3)
inner_outer_weights_summary_df[, 3] <- round(inner_outer_weights_summary_df[, 3], digits = 3)
inner_outer_weights_summary_df[, 4] <- as.character(percent(inner_outer_weights_summary_df[, 4], digits = 2))
inner_outer_weights_summary_df[, 5] <- as.character(percent(inner_outer_weights_summary_df[, 5], digits = 2))

# re-order the data frame
our_algorithms <- c("hungarian_algorithm", algorithm_names[grepl("blossom", algorithm_names)])
other_algorithms <- c("fast_greedy", "leading_eigen", "edge_betweenness",
                      "label_propagation", "louvain_algorithm", "walktrap_algorithm")
reordered_algorithm_names <- c(our_algorithms, other_algorithms)
reordered_inner_outer_weights_summary_df <- inner_outer_weights_summary_df[reordered_algorithm_names, ]

# output the data frame.
write.csv(reordered_inner_outer_weights_summary_df, file = output_path, quote = F, row.names = T)


