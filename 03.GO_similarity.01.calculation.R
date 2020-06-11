rm(list=ls())
ptm <- proc.time()
options(stringsAsFactors = F)
args<-commandArgs(TRUE)
library("magrittr")
library("org.Hs.eg.db")

#------------------source functions below---------------------------------------

conversion_table_path.TCGA <- file.path(getwd(), "function_scripts_directory", "mRNA_names_and_Entrez_IDs_dictionary.TCGA.csv")

source(file.path(getwd(), "function_scripts_directory", "GO_similarity_functions.R"))
#-------------------------------------------------------------------------------

if (is.na(args[1])) {
    # look up the input directory to get the cancer types.
    TCGA_correlation_coefficient_tables_path <- file.path(getwd(), "input_data_directory", "TCGA_correlation_coefficient_tables")
    
    cancer_types <- dir(TCGA_correlation_coefficient_tables_path) %>%
        strsplit(., "\\.") %>%
        sapply(., function(x) x[1][[1]])
} else {
    # use command line to run different cancer types in the backend concurrently.
    cancer_types <- args[1]
}


# Build the GO similarity databases.
ontology_terms <- c("BP", "CC", "MF")
similarity_database_list <- list()

for (the_ontology_term in ontology_terms) {
    sim_db <- godata('org.Hs.eg.db', ont=the_ontology_term, computeIC=TRUE)
    similarity_database_list[[the_ontology_term]] <- sim_db
}

for (cancer_type in cancer_types) {
    # cancer_type <- "BRCA"
    path_to_input_files <- file.path(getwd(),
                                     "output_data_directory",
                                     "binary_object_files",
                                     paste("communities_files", cancer_type, sep = "."))
    
    
    # get the names of different algorithms.
    communities_rounds <- dir(path_to_input_files) %>%
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
    
    
    # sapply(a_list_of_the_communities, modularity)
    
    output_file_names <- paste("GO_similarity", communities_rounds, "csv", sep = ".")
    path_to_output_files <- file.path(getwd(),
                                      "output_data_directory", 
                                      "plain_text_files",
                                      paste("GO_similarity_distance", cancer_type, sep = "."),
                                      output_file_names)
    
    dir.create(file.path(getwd(),
                         "output_data_directory",
                         "plain_text_files",
                         paste("GO_similarity_distance", cancer_type, sep = ".")),
               showWarnings = F, recursive = T)
    

    for (i in seq_along(communities_rounds)) {
        # i <- 12
        print(i)
        
        current_communities <- readRDS(path_to_input_files[i])
        cat("the communities produced from", communities_rounds[i], "has", length(communities(current_communities)), "clusters", sep = " ", "\n")

        # If the communities has only one cluster, it is meaningless to calculate the intra-cluster similarity scores or inter-cluster cluster similarity scores.
        if (length(communities(current_communities)) < 2) {
            cat(communities_rounds[i], "has only one cluster, so it is skipped.", sep = " ", "\n")
            next
        } else {
            GO_similarity_df <- calculate_GO_similiarity_from_a_communities(current_communities, conversion_table_path.TCGA, similarity_database_list)
            write.csv(GO_similarity_df, file = path_to_output_files[i], row.names = F, quote = F)
            cat(communities_rounds[i], "has its GO similarity scores calculated.", sep = " ", "\n")
        }
        # break
    }
    cat(cancer_type, "has done for GO similarity calculation.", sep = " ", "\n")
    print("******************************************************************")
    # break
}


print("Everything is done!")
proc.time() - ptm

