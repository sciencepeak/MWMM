library("igraph")
library("hash")
library("magrittr")
library("GOSemSim")
library("org.Hs.eg.db")


convert_a_communities_to_a_list_of_mRNA_only_clusters <- function(a_communites, conversion_table_path) {
    
    # now make a dictionary, using mRNA name as key and Entrez ID as value.
    name_ID_dictionary <- read.csv(conversion_table_path, stringsAsFactors = F) %$%
        hash(keys= mRNA_names, values = Entrez_IDs)
    
    get_Entrez_IDs_from_mRNA_names <- function(a_cluster_nodes) {
        a_cluster_mRNAs <- a_cluster_nodes[!grepl("hsa-", a_cluster_nodes)]
        a_cluster_Entrez_IDs <- hash::values(name_ID_dictionary, keys = a_cluster_mRNAs)
        return(a_cluster_Entrez_IDs)
    }
    a_list_of_clusters <- lapply(communities(a_communites), get_Entrez_IDs_from_mRNA_names)
    return(a_list_of_clusters)
}


calculate_average_GO_similarity_scores_in_one_cluster <- function(an_Entrez_ID_vector, a_sim_db, a_measure_method) {
    # an_Entrez_ID_vector <- c("25943", "80760", "79674")
    # clusterProfiler::bitr(an_Entrez_ID_vector, 'ENTREZID', 'GO', OrgDb='org.Hs.eg.db')
    # a_sim_db <- sim_db <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
    # a_measure_method <- "Resnik"
    
    similarity_matrix_trial <- try(mgeneSim(an_Entrez_ID_vector, semData = a_sim_db,
                                            measure = a_measure_method, drop = "IEA",
                                            combine = "avg", verbose = TRUE))
    
    # if there is no error.
    if (is.null(attr(similarity_matrix_trial, "condition")$message)) {
        similarity_matrix <- similarity_matrix_trial
        if (class(similarity_matrix) == "matrix" & length(similarity_matrix) == 0) {
            similarity_matrix <- 0
        }
        
    } else { # if there is error
        print("error found")
        similarity_matrix <- 0
    }
    
    if (length(an_Entrez_ID_vector) == 1) {
    	print("Only one gene in the cluster, the cluster is garbage.")
        similarity_matrix <- 0
    }
    
    similarity_score <- combineScores(similarity_matrix, combine = "avg")
    print(similarity_score)
    if (is.na(similarity_score)) {
        print(an_Entrez_ID_vector)
        stop("NA is found")
    }
    

    return(similarity_score)
}


calculate_average_GO_similarity_scores_among_clusters <- function(the_list_of_clusters, a_sim_db, a_measure_method) {
    # an_Entrez_ID_vector <- c("25943", "80760", "79674")
    # clusterProfiler::bitr(an_Entrez_ID_vector, 'ENTREZID', 'GO', OrgDb='org.Hs.eg.db')
    # a_sim_db <- sim_db <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
    # a_measure_method <- "Resnik"

    system.time(
        intercluster_score_matrix <- mclusterSim(the_list_of_clusters,
                                                 semData = a_sim_db,
                                                 measure = a_measure_method,
                                                 drop = "IEA",
                                                 combine = "avg")
    )

    inter_cluster_similarity_score <- combineScores(intercluster_score_matrix, combine = "avg")
    return(inter_cluster_similarity_score)
}


calculate_GO_similiarity_from_a_communities <- function(a_communities, path_to_the_conversion_table, a_GO_sim_db_list) {
    
    ontology_terms <- c("BP", "CC", "MF")
    measure_methods <- c("Resnik", "Lin", "Rel", "Jiang","Wang" )
    
    combination_df <- expand.grid(GO_terms = ontology_terms,
                                  measurements = measure_methods,
                                  stringsAsFactors = FALSE)
    ordered_combined_df <- combination_df[order(combination_df$GO_terms, combination_df$measurements), ]
    rownames(ordered_combined_df) <- seq_len(nrow(ordered_combined_df))

    the_list_of_clusters <- convert_a_communities_to_a_list_of_mRNA_only_clusters(a_communities, path_to_the_conversion_table)
    

    accum_vector_list <- list()
    for(a_row in seq_len(nrow(ordered_combined_df))){
        
        the_ontology_term <- ordered_combined_df[a_row, "GO_terms"]
        the_measure_method <- ordered_combined_df[a_row, "measurements"]
        sim_db <- a_GO_sim_db_list[[the_ontology_term]]
        
        GO_similarity_score_of_in_clusters <-
            sapply(the_list_of_clusters,
                   calculate_average_GO_similarity_scores_in_one_cluster,
                   a_sim_db = sim_db, a_measure_method = the_measure_method)
        
        # The average of the in-cluster similarity scores are weighted by cluster size.
        cluster_weights <- sapply(the_list_of_clusters, length)/sum(sapply(the_list_of_clusters, length))
        in_cluster_similarity_score <- weighted.mean(GO_similarity_score_of_in_clusters, cluster_weights)
        
        Sys.time()
        cat("intra-cluster scores of ", the_ontology_term, the_measure_method, "has been calculated.", sep = " ", "\n")

        inter_cluster_similarity_score <- calculate_average_GO_similarity_scores_among_clusters(the_list_of_clusters, sim_db, the_measure_method)
        
        Sys.time()
        cat("inter-cluster scores of ", the_ontology_term, the_measure_method, "has been calculated.", sep = " ", "\n")

        similarity_vector <- c(the_ontology_term,
                               the_measure_method,
                               round(in_cluster_similarity_score, digits = 3),
                               inter_cluster_similarity_score)
        names(similarity_vector) <- c("ontology",
                                      "measure",
                                      "in-cluster_similarity_score",
                                      "inter-cluster_similarity_score")
        
        accum_vector_list[[a_row]] <- similarity_vector
    }
    
    final_report_df <- as.data.frame(do.call(rbind, accum_vector_list))

    return(final_report_df)
}


#--------------------------------------Garbage code---------------------------
#------------------------------------non-parallel version total function.------------
# calculate_GO_similiarity_from_a_communities <- function(a_communities, path_to_the_conversion_table, a_GO_sim_db_list) {
#     
#     # a_communities <- current_communities
#     # path_to_the_conversion_table <- conversion_table_path.TCGA
#     # a_GO_sim_db_list <- similarity_database_list
# 
#     # -----------------------------------------------------------------------------
#     # --------------------- Set up parallel running below: ------------------------
#     library("foreach")
#     library("doParallel") # the parallel package will be loaded by doParallel
#     
#     # parameters to change on different system and machines:
#     # extra_libraries_path_directory <- "/net/home/lding/R/x86_64-redhat-linux-gnu-library/3.4/"
#     
#     # my packages are also installed on the following directories.
#     # .libPaths(c(.libPaths(), extra_libraries_path_directory))
#     
#     #setup parallel backend to use many processors, use all the available cores
#     available_cores_number <- detectCores()
#     if (available_cores_number > 15) {
#         available_cores_number <- 15
#     }
#     
#     parallel_cluster <- makeCluster(available_cores_number, outfile = file.path(getwd(), "stdout_std_err_for_parallel.txt"))
#     # session_information <- sessionInfo()
#     # if (session_information$platform == "x86_64-redhat-linux-gnu (64-bit)") {
#     #     parallel_cluster <- makeCluster(available_cores_number, type = "FORK")
#     # } else {
#     #     parallel_cluster <- makeCluster(available_cores_number, type = "PSOCK")
#     # }
#     # 
#     registerDoParallel(parallel_cluster)
#     
#     # load packages on each worker, so that all the threads can use packages by parLapply()
#     # http://stackoverflow.com/questions/17196261/understanding-the-differences-between-mclapply-and-parlapply-in-r
#     clusterEvalQ(parallel_cluster, library("igraph"))
#     clusterEvalQ(parallel_cluster, library("GOSemSim"))
#     clusterEvalQ(parallel_cluster, library("org.Hs.eg.db"))
#     
#     
#     clusterExport(cl = parallel_cluster, deparse(substitute(calculate_average_GO_similarity_scores_in_one_cluster)))
#     clusterExport(cl = parallel_cluster, deparse(substitute(calculate_average_GO_similarity_scores_among_clusters)))
#     
# 
#     ontology_terms <- c("BP", "CC", "MF")
#     measure_methods <- c("Resnik", "Lin", "Rel", "Jiang","Wang" )
#     
#     combination_df <- expand.grid(GO_terms = ontology_terms,
#                                   measurements = measure_methods,
#                                   stringsAsFactors = FALSE)
#     ordered_combined_df <- combination_df[order(combination_df$GO_terms, combination_df$measurements), ]
#     rownames(ordered_combined_df) <- seq_len(nrow(ordered_combined_df))
#     
#     clusterExport(cl = parallel_cluster, list("ordered_combined_df"), envir = environment())
#     # clusterExport(cl = parallel_cluster, list("ordered_combined_df", "a_GO_sim_db_list"), envir = environment())
# 
#     # final_report_df <- as.data.frame(t(result_df))
#     # return(final_report_df)
#     
#     
# 	the_list_of_clusters <- convert_a_communities_to_a_list_of_mRNA_only_clusters(a_communities, path_to_the_conversion_table)
# 
# 	# the_list_of_clusters <- the_list_of_clusters[1:5]
# 	# ontology_terms <- c("BP", "CC", "MF")[2]
# 	# measure_methods <- c("Resnik", "Lin", "Rel", "Jiang","Wang" )[2]
# 	
# 
# 	accum_vector_list <- foreach(a_row = seq_len(nrow(ordered_combined_df))) %dopar% {
# 
# 	    the_ontology_term <- ordered_combined_df[a_row, "GO_terms"]
# 	    the_measure_method <- ordered_combined_df[a_row, "measurements"]
# 	    sim_db <- a_GO_sim_db_list[[the_ontology_term]]
#         
# 	    GO_similarity_score_of_in_clusters <-
# 	        sapply(the_list_of_clusters,
# 	               calculate_average_GO_similarity_scores_in_one_cluster,
# 	               a_sim_db = sim_db, a_measure_method = the_measure_method)
# 	    
# 	    # The average of the in-cluster similarity scores are weighted by cluster size.
# 	    cluster_weights <- sapply(the_list_of_clusters, length)/sum(sapply(the_list_of_clusters, length))
# 	    in_cluster_similarity_score <- weighted.mean(GO_similarity_score_of_in_clusters, cluster_weights)
# 	    
# 	    Sys.time()
# 	    cat("intra-cluster scores of ", the_ontology_term, the_measure_method, "has been calculated.", sep = " ", "\n")
# 	    memory.size(max = TRUE)
# 	    gc()
# 	    # in_cluster_similarity_score <- 0.5
# 	    
# 	    inter_cluster_similarity_score <- calculate_average_GO_similarity_scores_among_clusters(the_list_of_clusters, sim_db, the_measure_method)
# 	    
# 	    Sys.time()
# 	    cat("inter-cluster scores of ", the_ontology_term, the_measure_method, "has been calculated.", sep = " ", "\n")
# 	    memory.size(max = TRUE)
# 	    gc()
# 	    
# 	    similarity_vector <- c(the_ontology_term,
# 	                           the_measure_method,
# 	                           round(in_cluster_similarity_score, digits = 3),
# 	                           inter_cluster_similarity_score)
# 	    names(similarity_vector) <- c("ontology",
# 	                                  "measure",
# 	                                  "in-cluster_similarity_score",
# 	                                  "inter-cluster_similarity_score")
# 	    
# 	    similarity_vector
# 	}
# 	
# 	final_report_df <- as.data.frame(do.call(rbind, accum_vector_list))
# 
# 	
# 	stopCluster(parallel_cluster)
# 	
# 	return(final_report_df)
# }




# 	
# 
#     ontology_terms <- c("BP", "CC", "MF")
#     measure_methods <- c("Resnik", "Lin", "Rel", "Jiang","Wang" )
#     accum_vector_list <- vector(mode = "list",
#                                 length = length(ontology_terms) * length(measure_methods))
# 
# 
#     tally_register <- 0
# 	for (the_ontology_term in ontology_terms) {
# 	    # construct the GO term database.
# 	    sim_db <- a_GO_sim_db_list[[the_ontology_term]]
# 	    for(the_measure_method in measure_methods) {
#             cat("The current GO term and mesurement method are ", the_ontology_term, the_measure_method, sep = " ", "\n")
# 	        system.time(
# 	        GO_similarity_score_of_in_clusters <-
# 	            sapply(the_list_of_clusters,
# 	                   calculate_average_GO_similarity_scores_in_one_cluster,
# 	                   a_sim_db = sim_db, a_measure_method = the_measure_method)
# 	        )
# 
# 	        # The average of the in-cluster similarity scores are weighted by cluster size.
# 	        cluster_weights <- sapply(the_list_of_clusters, length)/sum(sapply(the_list_of_clusters, length))
# 	        in_cluster_similarity_score <- weighted.mean(GO_similarity_score_of_in_clusters, cluster_weights)
#             
# 	        # in_cluster_similarity_score <- 0.5
# 	        
# 	        inter_cluster_similarity_score <- calculate_average_GO_similarity_scores_among_clusters(the_list_of_clusters, sim_db, the_measure_method)
# 
# 	        
# 
# 
# 	        similarity_vector <- c(the_ontology_term,
# 	                           the_measure_method,
# 	                           round(in_cluster_similarity_score, digits = 3),
# 	                           inter_cluster_similarity_score)
# 	        names(similarity_vector) <- c("ontology",
# 	                                  "measure",
# 	                                  "in-cluster_similarity_score",
# 	                                  "inter-cluster_similarity_score")
# 
# 	        tally_register <- tally_register + 1
# 	        cat("now", tally_register, "has been finished", "\n", sep = " ")
# 	        accum_vector_list[[tally_register]] <- similarity_vector
# 	    }
# 	}
# 
# 	final_report_df <- as.data.frame(do.call(rbind, accum_vector_list))
# 	
# 	stopCluster(parallel_cluster)
# 	
# 	return(final_report_df)
# }

# ----------------------------test-----------------------------

# parallel_GO_function <- function(a_row, the_list_of_clusters) {
#     
#     the_ontology_term <- a_row[1]
#     the_measure_method <- a_row[2]
#     sim_db <- godata('org.Hs.eg.db', ont=the_ontology_term, computeIC=TRUE)
#     
#     system.time(
#         GO_similarity_score_of_in_clusters <-
#             sapply(the_list_of_clusters,
#                    calculate_average_GO_similarity_scores_in_one_cluster,
#                    a_sim_db = sim_db, a_measure_method = the_measure_method)
#     )
#     
#     # The average of the in-cluster similarity scores are weighted by cluster size.
#     cluster_weights <- sapply(the_list_of_clusters, length)/sum(sapply(the_list_of_clusters, length))
#     in_cluster_similarity_score <- weighted.mean(GO_similarity_score_of_in_clusters, cluster_weights)
# 
#     
#     # system.time(
#     #     among_cluster_similarity_matrix <- mclusterSim(the_list_of_clusters,
#     #                                                    semData = sim_db,
#     #                                                    measure = the_measure_method,
#     #                                                    drop = "IEA", combine = "avg")
#     # )
#     # inter_cluster_similarity_score <- combineScores(among_cluster_similarity_matrix, combine = "avg")
#     # inter_cluster_similarity_score
#     
#     inter_cluster_similarity_score <- calculate_average_GO_similarity_scores_among_clusters(the_list_of_clusters, sim_db, the_measure_method)
# 
#     
#     similarity_vector <- c(the_ontology_term,
#                            the_measure_method,
#                            round(in_cluster_similarity_score, digits = 3),
#                            inter_cluster_similarity_score)
#     names(similarity_vector) <- c("ontology",
#                                   "measure",
#                                   "in-cluster_similarity_score",
#                                   "inter-cluster_similarity_score")
#     return(similarity_vector)
# }
# 
# 
# calculate_GO_similiarity_from_a_communities2 <- function(a_communities, path_to_the_conversion_table) {
#     
#     
#     # -----------------------------------------------------------------------------
#     # --------------------- Set up parallel running below: ------------------------
#     library("foreach")
#     library("doParallel") # the parallel package will be loaded by doParallel
#     
#     # parameters to change on different system and machines:
#     # extra_libraries_path_directory <- "/net/home/lding/R/x86_64-redhat-linux-gnu-library/3.4/"
#     
#     # my packages are also installed on the following directories.
#     # .libPaths(c(.libPaths(), extra_libraries_path_directory))
#     
#     #setup parallel backend to use many processors, use all the available cores
#     available_cores_number <- detectCores()
#     if (available_cores_number > 15) {
#         available_cores_number <- 15
#     }
#     session_information <- sessionInfo()
#     if (session_information$platform == "x86_64-redhat-linux-gnu (64-bit)") {
#         parallel_cluster <- makeCluster(available_cores_number, type = "FORK")
#     } else {
#         parallel_cluster <- makeCluster(available_cores_number, type = "PSOCK")
#     }
#     
#     registerDoParallel(parallel_cluster)
#     
#     # load packages on each worker, so that all the threads can use packages by parLapply()
#     # http://stackoverflow.com/questions/17196261/understanding-the-differences-between-mclapply-and-parlapply-in-r
#     clusterEvalQ(parallel_cluster, library("igraph"))
#     clusterEvalQ(parallel_cluster, library("GOSemSim"))
#     clusterEvalQ(parallel_cluster, library("org.Hs.eg.db"))
#     
#     
#     clusterExport(cl = parallel_cluster, deparse(substitute(calculate_average_GO_similarity_scores_in_one_cluster)))
#     clusterExport(cl = parallel_cluster, deparse(substitute(calculate_average_GO_similarity_scores_among_clusters)))
#     
#     current_list_of_clusters <- convert_a_communities_to_a_list_of_mRNA_only_clusters(a_communities, path_to_the_conversion_table)
#     # load("a_communities_of_star_graphs.markov_matrix.RData")
#     # the_list_of_clusters <- convert_a_communities_to_a_list_of_mRNA_only_clusters(a_communities_of_star_graphs.markov_matrix, "mRNA_names_and_Entrez_IDs_dictionary.csv")
#     
#     # current_list_of_clusters <- current_list_of_clusters[1:2]
#     
#     ontology_terms <- c("BP", "CC", "MF")
#     measure_methods <- c("Resnik", "Lin", "Rel", "Jiang","Wang" )
#     
#     combination_df <- expand.grid(GO_terms = ontology_terms,
#                                   measurements = measure_methods)
#     ordered_combined_df <- combination_df[order(combination_df$GO_terms, combination_df$measurements), ]
#     result_df <- parApply(cl = parallel_cluster, ordered_combined_df, 1, parallel_GO_function, the_list_of_clusters = current_list_of_clusters)
#     
#     
#     stopCluster(parallel_cluster)
#     
#     final_report_df <- as.data.frame(t(result_df))
#     return(final_report_df)
# }

# calculate_average_GO_similarity_scores_among_clusters <- function(the_list_of_clusters, a_sim_db, a_measure_method) {
#     # an_Entrez_ID_vector <- c("25943", "80760", "79674")
#     # clusterProfiler::bitr(an_Entrez_ID_vector, 'ENTREZID', 'GO', OrgDb='org.Hs.eg.db')
#     # a_sim_db <- sim_db <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
#     # a_measure_method <- "Resnik"
#     
#     system.time(
#         intercluster_score_matrix <- mclusterSim(the_list_of_clusters,
#                                                  semData = a_sim_db,
#                                                  measure = a_measure_method,
#                                                  drop = "IEA",
#                                                  combine = "avg")
#     )
#     
#     # cluster_size_vector <- sapply(the_list_of_clusters, length)
#     # remaining_cluster_names <- colnames(intercluster_score_matrix)
#     # remaining_cluster_sizes <- cluster_size_vector[remaining_cluster_names]
#     # remaining_cluster_number <- length(remaining_cluster_names)
#     # 
#     # 
#     # cluster_size_matrix <- matrix(data = 0,
#     #                               nrow = remaining_cluster_number,
#     #                               ncol = remaining_cluster_number,
#     #                               dimnames = list(rownames(intercluster_score_matrix),
#     #                                               colnames(intercluster_score_matrix)))
#     # 
#     # for (i in rownames(cluster_size_matrix)) {
#     #     for (j in colnames(cluster_size_matrix)) {
#     #         cluster_size_matrix[i, j] <- cluster_size_vector[as.numeric(i)] + cluster_size_vector[as.numeric(j)]
#     #     }
#     # }
#     # 
#     # weighted_size_matrix <- cluster_size_matrix / 2 / sum(remaining_cluster_sizes) / remaining_cluster_number
#     # # sum(weighted_size_matrix)
#     # 
#     # cat("the length of intercluster score matrix is", length(intercluster_score_matrix), sep = " ", "\n")
#     # # print(intercluster_score_matrix)
#     # cat("the length of weighted size matrix is ", length(weighted_size_matrix), sep = " ", "\n")
#     # # print(weighted_size_matrix)
#     # # print(the_list_of_clusters)
#     # print("---------------------------------------------------------")
#     # # sum(weighted_size_matrix)
#     # # final_matrix <- intercluster_score_matrix * weighted_size_matrix
#     # # sum(final_matrix)
#     # inter_cluster_similarity_score <- weighted.mean(intercluster_score_matrix, weighted_size_matrix)
#     # 
#     # inter_cluster_similarity_score <- round(inter_cluster_similarity_score, digits = 3)
#     
#     
#     inter_cluster_similarity_score <- combineScores(intercluster_score_matrix, combine = "avg")
#     
#     
#     return(inter_cluster_similarity_score)
# }
