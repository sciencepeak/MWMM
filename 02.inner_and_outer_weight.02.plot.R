rm(list = ls())
ptm <- proc.time()
options(stringsAsFactors = F)
library("magrittr")

# look up the input directory to get the cancer types.
TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])

graph_output_directory <- file.path(getwd(),
                                    "output_data_directory",
                                    "graph_files",
                                    "inner_and_outer_weight_plots")
dir.create(graph_output_directory, showWarnings = F, recursive = T)



only_draw_MWMM <- TRUE
# only_draw_MWMM <- FALSE

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
    other_communities_rounds <- c("label_propagation",
                                  "edge_betweenness",
                                  "walktrap_algorithm",
                                  "fast_greedy",
                                  "leading_eigen",
                                  "louvain_algorithm")
    
    if (only_draw_MWMM == TRUE) {
        communities_rounds <- our_communities_rounds
    } else {
        # draw graph of all algorithm results.
        communities_rounds <- c(our_communities_rounds, other_communities_rounds)
    }
    
    input_file_names <- paste("inner_and_outer_weight", communities_rounds, "csv", sep = ".")
    path_to_input_files <- file.path(getwd(),
                                     "output_data_directory",
                                     "plain_text_files",
                                     paste("inner_and_outer_weight", cancer_type, sep = "."),
                                     input_file_names)
    
    
    a_list_of_algorithm_condition_weights <- lapply(path_to_input_files, read.csv, check.names = FALSE)
    names(a_list_of_algorithm_condition_weights) <- communities_rounds
    
    data_table <- data.frame(algorithm_round = character(0),
                             cluster_number = numeric(0),
                             condition_1_Yes_percent = numeric(0),
                             condition_2_Yes_percent = numeric(0))
    
    get_true_percent <- function(a_column_vector) {
        x_percent <- sum(a_column_vector) / length(a_column_vector)
        rounded_percent <- round(x_percent, digits = 2)
        rounded_percent
    }
    
    for (i in seq_along(a_list_of_algorithm_condition_weights)) {
        a_condition_df <- a_list_of_algorithm_condition_weights[[i]]
        
        data_table[i, "algorithm_round"] <- communities_rounds[i]
        data_table[i, "cluster_number"] <- nrow(a_condition_df)
        data_table[i, "condition_1_Yes_percent"] <- get_true_percent(a_condition_df$condition_01)
        data_table[i, "condition_2_Yes_percent"] <- get_true_percent(a_condition_df$condition_02)
        
    }
    
    
    # list our algorithm first and based on merging order.
    # list other algorithms following and based on the decrease of the cluster number.
    our_algorithm_df <- data_table[!data_table$algorithm_round %in% other_communities_rounds, ]
    other_algorithm_df <- data_table[data_table$algorithm_round %in% other_communities_rounds, ]
    
    sorted_other_algorithm_df <- other_algorithm_df[order(-other_algorithm_df$cluster_number), ]
    data_table <-rbind.data.frame(our_algorithm_df, sorted_other_algorithm_df)
    
    # data_table <- data_table[order(-data_table$cluster_number), ]
    
    algorithm_number <- length(communities_rounds)
    
    
    
    if (only_draw_MWMM == TRUE) {
        path_to_output_file <- file.path(graph_output_directory, paste("inner_and_outer_weight_conditions.MWMM_only", cancer_type, "pdf", sep = "."))
    } else {
        path_to_output_file <- file.path(graph_output_directory, paste("inner_and_outer_weight_conditions", cancer_type, "pdf", sep = "."))
    }

    pdf(path_to_output_file)
    
    my_colors <- c("red", "green", "steelblue")
    my_pchs <- c(15, 16, 17)
    
    # mar:numerical vector indicating margin size c(bottom, left, top, right) in lines. default = c(5, 4, 4, 2) + 0.1 
    
    par(pin = c(6.69, 6.69), mar=c(10, 5, 5, 5), xpd=TRUE)
    
    plot(c(1, algorithm_number), c(1, 350), type = "n", frame.plot = F,
         xaxt="n", xlab = "", ylab = "cluster number (count)")
    
    axis(side = 1, at=1:algorithm_number, labels= data_table$algorithm_round, las = 2)
    lines(1:algorithm_number, data_table$cluster_number,
          type = "o", lwd = 2, cex = 1.5,
          col = my_colors[1], pch = my_pchs[1])
    
    mtext("algorithm name", side = 1, line = 8)
    
    par(new = T)
    cond_1 <- data_table$condition_1_Yes_percent * 100
    cond_2 <- data_table$condition_2_Yes_percent * 100
    
    plot(1:algorithm_number, cond_1, axes = F, xlab = NA, ylab = NA, ylim = c(1, 100),
         type = "o", lwd = 2, cex = 1.5,
         col = my_colors[2], pch = my_pchs[2])
    
    par(new = T)
    plot(1:algorithm_number, cond_2, axes = F, xlab = NA, ylab = NA, ylim = c(1, 100),
         type = "o", lwd = 2, cex = 1.5,
         col = my_colors[3], pch = my_pchs[3])
    
    right_y_axis_ticks <- seq(0, 100, by = 10)
    axis(side = 4, at=right_y_axis_ticks,  labels = right_y_axis_ticks)
    mtext("meeting condition (%)", side = 4, line = 3)
    
    
    legend("topleft", inset = c(0.01, -0.23),
           legend = c("cluster number",
                      "condition 1",
                      "condition 2"),
           col = my_colors, pch = my_pchs, lty = 1, lwd = 2)
    
    dev.off()

}


proc.time() - ptm
