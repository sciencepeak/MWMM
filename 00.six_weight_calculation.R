rm(list = ls())
ptm <- proc.time()
options(stringsAsFactors = F)
library("magrittr")
library("hash")

# gene symbol correction. The dictionary table is adapted from package "HGNChelper"
# build a dictionary for mogrified query date-style gene symbol.
mo_map_df <- read.csv(file.path(getwd(), "function_scripts_directory", "mog_map.csv"), check.names = F)
mo_map_dictionary <- hash(keys = mo_map_df$mogrified, values = mo_map_df$original)

# create the path for the output data in a platform-independent way.

dir.create(file.path(getwd(), "output_data_directory", "binary_object_files"),
           showWarnings = F, recursive = T)
dir.create(file.path(getwd(), "output_data_directory", "graph_files", "six_weights_distribution"),
           showWarnings = F, recursive = T)
dir.create(file.path(getwd(), "output_data_directory", "plain_text_files"),
           showWarnings = F, recursive = T)


# look up the input directory to get the cancer types.
TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")

cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
    strsplit(., "\\.") %>%
    sapply(., function(x) x[1][[1]])


# read the cancer correlation data file and correct Excel-date corrupted gene symbols
a_list_of_input_df <- list()

for (cancer_type in cancer_types) {

    # create the output directory for calculated weights edge list.
    dir.create(file.path(getwd(), "output_data_directory",
                         "plain_text_files",
                         paste("six_kinds_of_weights_files", cancer_type, sep = ".")),
               showWarnings = F, recursive = T)
    
    # read the correlation data file into a data frame.
    input_file_name <- file.path(TCGA_correlation_coefficient_tables_path,
                                 paste(cancer_type, "csv", sep = "."))
    input_df <- read.csv(input_file_name)
    
    temporary_gene_names <- input_df$mRNA
    
    index_mogrified <- temporary_gene_names %in% mo_map_df$mogrified
    
    # if there are match of the Excel date mogrified gene symbols.
    if (sum(index_mogrified) > 0) {
        
        corrected_gene_names <- hash::values(mo_map_dictionary, keys = temporary_gene_names[index_mogrified])
        
        temporary_gene_names[index_mogrified] <- corrected_gene_names
        
        input_df$mRNA <- temporary_gene_names
        
        a_list_of_input_df[[cancer_type]] <- input_df
        
    } else { # if there is no weired gene symbols.
        a_list_of_input_df[[cancer_type]] <- input_df
    }
}



# --------------------------------------------------------------------------
# Calculate the six kinds of proposed weights for all cancer types.
for (cancer_type in names(a_list_of_input_df)) {
    input_df <- a_list_of_input_df[[cancer_type]]
    
    # The red edges denote that the correlation coefficients are converted from negative in normal to positive in tumor.
    # The green edges denote that the correlation coefficients are converted from positive in normal to negative in tumor.
    correlation_change <- rep("green", times = nrow(input_df))
    correlation_change[input_df$T_CC > 0] <- "red"
    
    # get different weight calculated.
    
    # either T_CC or N_CC, as long as it's negative, we take it as the weight
    all_negative_value_weight <- numeric(nrow(input_df))
    all_negative_value_weight[input_df$T_CC < 0] <- input_df$T_CC[input_df$T_CC < 0]
    all_negative_value_weight[input_df$N_CC < 0] <- input_df$N_CC[input_df$N_CC < 0]
    all_negative_value_weight <- abs(all_negative_value_weight)
    
    # either T_CC or N_CC, as long as it's positive, we take it as the weight
    all_positive_value_weight <- numeric(nrow(input_df))
    all_positive_value_weight[input_df$T_CC >0] <- input_df$T_CC[input_df$T_CC > 0]
    all_positive_value_weight[input_df$N_CC > 0] <- input_df$N_CC[input_df$N_CC > 0]
    
    # the sum of absolute value of T_CC and N_CC is the weight
    arithmetic_mean_value_weight <- (abs(input_df$T_CC) + abs(input_df$N_CC))/2
    
    # the product of absolute value of T_CC and N_CC is the weight
    geometric_mean_value_weight <- sqrt(abs(input_df$T_CC * input_df$N_CC))
    
    
    # integrated mean needs to calculate different lambda
    
    mean_T_CC_positive <- mean(input_df$T_CC[input_df$T_CC > 0])
    mean_T_CC_negative <- abs(mean(input_df$T_CC[input_df$T_CC < 0]))
    mean_N_CC_positive <- mean(input_df$N_CC[input_df$N_CC > 0])
    mean_N_CC_negative <- abs(mean(input_df$N_CC[input_df$N_CC < 0]))
    
    lambda_1 <- mean_T_CC_positive/(mean_T_CC_positive + mean_N_CC_negative)
    one_minus_lambda_1 <- mean_N_CC_negative/(mean_T_CC_positive + mean_N_CC_negative)
    lambda_2 <- mean_T_CC_negative/(mean_T_CC_negative + mean_N_CC_positive)
    one_minus_lambda_2 <- mean_N_CC_positive/(mean_T_CC_negative + mean_N_CC_positive)
    
    integrated_mean_value_weight <- numeric(nrow(input_df))
    integrated_mean_value_weight[input_df$T_CC >0] <-
        lambda_1 * input_df$T_CC[input_df$T_CC > 0] + one_minus_lambda_1 * abs(input_df$N_CC[input_df$T_CC > 0])
    integrated_mean_value_weight[input_df$T_CC <0] <-
        lambda_2 * abs(input_df$T_CC[input_df$T_CC < 0]) + one_minus_lambda_2 * input_df$N_CC[input_df$T_CC < 0]
    
    # calculate the pairwise maximum between two numeric vectors
    maximum_absolute_value_weight <- pmax(abs(input_df$T_CC), abs(input_df$N_CC))
    
    
    # put all the weight vectors into a list to loop on them later
    all_weight_list <- list(integrated_mean_value_weight = integrated_mean_value_weight,
                            all_negative_value_weight = all_negative_value_weight,
                            all_positive_value_weight = all_positive_value_weight,
                            arithmetic_mean_value_weight = arithmetic_mean_value_weight,
                            geometric_mean_value_weight = geometric_mean_value_weight,
                            maximum_absolute_value_weight = maximum_absolute_value_weight)
    
    # assemble the edgelist and write them to the hard disk.
    for (a_weight_name in names(all_weight_list)) {
        
        current_weight_vector <- all_weight_list[[a_weight_name]]

        current_edgelist_dataframe <- data.frame(mRNA = input_df$mRNA,
                                                 microRNA = input_df$microRNA,
                                                 weight = current_weight_vector,
                                                 edge_color = correlation_change)
        current_output_file_name <- paste(a_weight_name, "csv", sep = ".")
        
        output_file <- file.path(getwd(), "output_data_directory", "plain_text_files",
                                 paste("six_kinds_of_weights_files", cancer_type, sep = "."),
                                 current_output_file_name)
        
        write.csv(current_edgelist_dataframe,
                  file = output_file,
                  quote = F, row.names = F)
    }

    
    # pause the execution for 1 second such that the later file can be written in a later time.
    Sys.sleep(1)
    

    pdf(file = file.path(getwd(), "output_data_directory", "graph_files", "six_weights_distribution",
                         paste("distribution_of_six_weights", cancer_type, "pdf", sep = ".")),
        width = 14, height = 21, pointsize = 14)
    
    # plot 3x2 graphs in one page.
    par(mfrow=c(3,2),
        mar=c(5.1, 5.1, 4.1, 2.1),
        mgp=c(3, 0.8, 0),
        cex.axis = 2,
        cex.lab = 2,
        cex.main = 2)
    
    # The for loop draws histogram graphs to show the distribution of the different weight.
    for (i in seq_along(all_weight_list)) {
        
        # the calcuated weight is so small like 0.3312909 that we need to round it.
        round_digits <- 3
        
        current_weight_vector <- all_weight_list[[i]]
        current_weight_name <- names(all_weight_list)[i]
        current_figure_title <- paste("Distribution of", current_weight_name,
                                      "rounded to", round_digits, "decimal places",
                                      sep = " ")
        current_figure_title <- paste("distribution of", current_weight_name, sep = " ")
        # To plot the histogram graph with proper breaks, we set up breaks
        # based on the length of table(current_weight_vector).
        # If the length of the frequency table is 350, we round it to 400 
        actual_breaks <- current_weight_vector %>%
            round(digits = round_digits) %>%
            table %>%
            length
        
        apt_breaks <- actual_breaks %>%
            round(digits = -2)
        
        print(actual_breaks)
        print(apt_breaks)
        
        # if the script rounds the actual breaks 740 to 700,
        # then we add 100 to 700 to make the apt breaks to be 800
        if (apt_breaks < actual_breaks) {
            apt_breaks <- apt_breaks + 100
        }
        print(apt_breaks)
        print("------------------------------------")
        
        # calculate the histogram
        frequency_histogram <- current_weight_vector %>%
            round(digits = round_digits) %>%
            hist(breaks = apt_breaks, col="black", xlab = "weight", ylab = "count",
                 main = current_figure_title)
        
        # # plot the histogram to the hard disk.
        # output_graph_name <- paste(current_figure_title, "png", sep = ".")
        # 
        # png(filename = file.path(result_files_directory, output_graph_name),
        #     width = 1200, height = 1200, bg = "transparent")
        # plot(frequency_histogram, col="red", xlab = "weight", ylab = "count",
        #      main = current_figure_title)
        # dev.off()
    }
    dev.off()
}

proc.time() - ptm
