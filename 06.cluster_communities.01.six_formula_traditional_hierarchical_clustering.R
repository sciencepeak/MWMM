rm(list = ls())
ptm <- proc.time()
options(stringsAsFactors = F)

# This script draw the traditional hierarchical clustering results of six weight formulas
# at certain steps n, here we choose n = 38, because that is the threshold of integrated mean value weights.

library("igraph")

source(file.path(getwd(), "function_scripts_directory", "draw_cluster.R"))
cancer_type <- "BRCA"

# number of steps the algorithm loop runs is based on the weight threshold.
weight_threshold = "desired"

# traditional hierarchical clustering method applies to different 6 formulas

if (weight_threshold == 0.5) steps = 378
if (weight_threshold == 0.6) steps = 85
if (weight_threshold == 0.7) steps = 38

if (weight_threshold == "desired") steps = 38


weights_directory <- file.path(getwd(),
                               "output_data_directory",
                               "plain_text_files",
                               paste("six_kinds_of_weights_files", cancer_type, sep = "."))


weights_graph_dir <- file.path(getwd(), "output_data_directory", "graph_files",
                               paste("six_weights_traditional_clustering", cancer_type, sep = "."))
dir.create(weights_graph_dir, showWarnings = F, recursive = T)
input_file_name_vector <- c("all_negative_value_weight.csv",
                            "all_positive_value_weight.csv",
                            "arithmetic_mean_value_weight.csv",
                            "geometric_mean_value_weight.csv",
                            "integrated_mean_value_weight.csv",
                            "maximum_absolute_value_weight.csv")

# Read the six formula weight files into a list.
all_formula_edgelist_list <- list()

for ( current_file_name in input_file_name_vector) {

    # the base name doesn't have the .csv suffix
    current_base_name <- gsub("\\.csv", "", current_file_name)
    current_path_to_file <- file.path(weights_directory, current_file_name)
    
    # input dataframe is a csv edgelist file with columns: microRNA, mRNA, weight, role
    input_df <- read.csv(current_path_to_file, check.names=FALSE)
    
    all_formula_edgelist_list[[current_base_name]] <- input_df
}

# sort and select the top weighted edges in different 6 formulas.

top_weight_six_formula_list <- list()
for (current_edgelist_name in names(all_formula_edgelist_list)) {
    
    current_edgelist_df <- all_formula_edgelist_list[[current_edgelist_name]]
    
    # sort the edgelist dataframe descreasingly by weight:  (highest -> lowest)
    sorted_edgelist_df <- current_edgelist_df[(order(current_edgelist_df$weight, decreasing = T)), ]
    
    # Only select the top n weights
    filtered_edgelist_df <- sorted_edgelist_df[1:steps, ]
    
    top_weight_six_formula_list[[current_edgelist_name]] <- filtered_edgelist_df
}

for (formula_name in names(all_formula_edgelist_list)) {
    one_raw_edgelist_df <- all_formula_edgelist_list[[formula_name]]
    one_formula_edgelist_df <- top_weight_six_formula_list[[formula_name]]
    explanation_part <- paste("top", steps, "edges", sep = "_")
    current_output_file <- paste(formula_name, explanation_part,cancer_type, "pdf", sep = ".")
    print(current_output_file)
    output_file <- file.path(weights_graph_dir, current_output_file)
    
    draw_cluster(one_formula_edgelist_df, one_formula_edgelist_df, output_file)
}


proc.time() - ptm
