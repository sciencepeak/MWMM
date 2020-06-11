options(stringsAsFactors = F)
library("igraph")
library("clue")
library("magrittr") # R pipeline operations.
library("reshape2")


calculate_hungarian_algorithm <- function(input_file_path) {

    raw_edge_list_df <- read.csv(input_file_path, check.names = F)[, 1:3]
    raw_edge_list_df <- raw_edge_list_df[raw_edge_list_df$mRNA != "SLC35E2", ]

    # ensure that there are no 0 in the weight.
    # if there are 0 in the weight, quit the whole program.
    if ( 0 %in% raw_edge_list_df$weight){
        print("There are zero weight.")
        quit(status = -1)
    }
    
    raw_graph <- graph.data.frame(raw_edge_list_df, directed = F)
    
    raw_adjacency_matrix <- as_adjacency_matrix(raw_graph, type="both",
                                               names=TRUE, sparse=FALSE, attr="weight")

    # The mRNA names is much more than microRNA in number.
    # however, solve_LSAP(raw_bipartite_graph_matrix, maximum = T) :
    # x must not have more rows than columns.
    # Thereby, the row names have to be microRNA, the column names have to be mRNA.
    

    raw_bipartite_graph_matrix <- raw_edge_list_df %>%
        # convert the edge list to the matrix, ie. long format to wide format.
        dcast(., microRNA ~ mRNA, value.var = "weight", fill = 0, fun.aggregate = mean) %$%
        
        # change the rownames to the first column values
        set_rownames(., microRNA) %>%
        
        # Remove the first column, since the first column values becomes row names.
        inset(., "microRNA", value = NULL) %>%
        
        # Convert from data frame to matrix type
        as.matrix

    #------the above block generates reference matrix-----------------
    
    # ---------------------Apply the Hungarian algorithms-------------------------
    # ensure that the number of rows is less than the number of columns.
    # this is required by the clue package::solve_LSAP
    if ( nrow(raw_bipartite_graph_matrix) >= ncol(raw_bipartite_graph_matrix)){
        print("number of rows >= numbe of columns in raw bipartite gaph matrix")
        quit(status = -1)
    }
    
    # Get the edge list of the mRNA leave nodes allocated to different miRNA center nodes.
    star_graphs_weighed_edge_list_df <- data.frame(microRNA = character(0), mRNA = character(0), weight = numeric(0))
    
    
    # intialize the current reference matrix that will be deducted columns.
    current_raw_bipartite_graph_matrix <- raw_bipartite_graph_matrix
    
    
    i = 1
    
    while (TRUE) {
        
        # intialize the output edge list data frame.
        edge_list_df <- data.frame(mRNA = character(0), microRNA = character(0), weight = numeric(0))
        
        # hungarian alogrithm on a matrix
        # refer to https://cran.r-project.org/web/packages/clue/clue.pdf
        
        current_hungarian_matching <- solve_LSAP(current_raw_bipartite_graph_matrix, maximum = T)
        
        # get the current hungarian matching vertexes and edges, dump them into an edge list data frame.
        for (j in seq_along(current_hungarian_matching)) {
            miRNA_names <- rownames(current_raw_bipartite_graph_matrix)[j]
            mRNA_names <- colnames(current_raw_bipartite_graph_matrix)[current_hungarian_matching[j]]
            matrix_entry <- current_raw_bipartite_graph_matrix[miRNA_names, mRNA_names]
            
            edge_list_df[j, "mRNA"] <- mRNA_names
            edge_list_df[j, "microRNA"] <- miRNA_names
            edge_list_df[j, "weight"] <- matrix_entry
        }
        
        # # only retain the non-zero weight matching.
        # missing_miRNAs %in% edge_list_df$microRNA
        
        # edge_list_df[edge_list_df$microRNA %in% missing_miRNAs, ]
        
        cat("before remove zero hungarian match, the edge list length is", nrow(edge_list_df), sep = " ", "\n")
        edge_list_df <- edge_list_df[edge_list_df$weight != 0, ]
        
        cat("current Hungarian matching get an edge list of",  nrow(edge_list_df), "rows.", sep = " ", "\n")
        # print(edge_list_df)
        # # missing_miRNAs %in% edge_list_df$microRNA
        # print("what is the bipartite graph before removing zero rows?")
        # print(current_raw_bipartite_graph_matrix)
    
        # update the current reference matrix that has columns deducted.
        remaining_column_names <- setdiff(colnames(current_raw_bipartite_graph_matrix), edge_list_df$mRNA)
        # cat("remaining mRNAs length is", length(remaining_column_names), sep = " ", "\n")

        current_raw_bipartite_graph_matrix <- current_raw_bipartite_graph_matrix[, remaining_column_names, drop=F]
        
        
        # missing_miRNAs %in% rownames(current_raw_bipartite_graph_matrix)
        # dim(current_raw_bipartite_graph_matrix)
        # rowSums(current_raw_bipartite_graph_matrix)
        
        row_index <- rowSums(current_raw_bipartite_graph_matrix) != 0
                current_raw_bipartite_graph_matrix <- current_raw_bipartite_graph_matrix[row_index, , drop=F]

        cat("the current bipartite graph matrix has dimension:", dim(current_raw_bipartite_graph_matrix), sep = " ", "\n")

        star_graphs_weighed_edge_list_df <- rbind(star_graphs_weighed_edge_list_df, edge_list_df)
        cat("the Hungarian algorithm has run", i, "rounds.", sep = " ", "\n")
        i = i + 1
        print("=================================================")
        if (length(current_raw_bipartite_graph_matrix) == 0) break
    }
    
    # Because some miRNAs might be removed from the final edge list, a new graph is needed.
    secondary_graph <- graph.data.frame(star_graphs_weighed_edge_list_df, directed = F)
    secondary_adjacency_matrix = as_adjacency_matrix(secondary_graph, type="both",
                                                     names=TRUE, sparse=FALSE, attr="weight")
    
    # get a list of star graphs represented by miRNA names and mRNA names.
    star_graph_list <- split(star_graphs_weighed_edge_list_df$mRNA,
                             list(star_graphs_weighed_edge_list_df$microRNA))
    
    # -------------------------------------------------------------------------
    # The old way of organizing the star graph objects is to put them into a list.
    
    a_list_of_star_graph_objects.hungarian_algorithm <- vector(mode = "list",
                                                               length = length(star_graph_list))
    
    # i from 1 to 312 in the case of BRCA.
    for(i in seq_along(star_graph_list)) {
        current_center_name <- names(star_graph_list)[i]
        current_leaf_names <- star_graph_list[[i]]
        
        # ensure that the star has at least one leaf.
        if ( length(current_leaf_names) == 0){
            quit(status = -1)
        }
        
        # make the star graph for the current row of the edge list data frame.
        current_star_graph_object <- make_star(length(current_leaf_names) + 1,
                                               mode = "out", center = 1)
        
        # assign the microRNA and mRNA names to the vertexes.
        V(current_star_graph_object)$name <- c(current_center_name, current_leaf_names)
        
        # append the star graph from one row to the star graph list for one edge list data frame.
        a_list_of_star_graph_objects.hungarian_algorithm[[i]] <- current_star_graph_object
    }
    
    names(a_list_of_star_graph_objects.hungarian_algorithm) <- names(star_graph_list)

    # ----------------------------------------------------------------------
    # The new way of organizing clusters is to construct a communities object.
    # The hardest job is to create a named vector as the membership of a communities object.
    
    all_miRNA_members <- seq_along(star_graph_list)
    names(all_miRNA_members) <- names(star_graph_list)
    
    all_mRNA_members <- numeric(0)
    for (i in seq_along(star_graph_list)) {
        current_center_name <- names(star_graph_list)[i]
        current_leaf_names <- star_graph_list[[i]]
        
        # ensure that the star has at least one leaf.
        if ( length(current_leaf_names) == 0){
            quit(status = -1)
        }
        current_mRNA_members <- rep(i, times = length(current_leaf_names))
        names(current_mRNA_members) <- current_leaf_names
        all_mRNA_members <- append(all_mRNA_members, current_mRNA_members)
    }
    
    # all_node_members is used for membership of the communities object.
    all_node_members <- append(all_mRNA_members, all_miRNA_members)
    cat("the number of nodes is: ", length(all_node_members), sep=" ", "\n")
    

    # order the membership based on the row names order of the raw adjacency matrix.
    membership_vector <- all_node_members[rownames(secondary_adjacency_matrix)]
    cat("the number of members is: ", length(membership_vector), sep=" ", "\n")
    
    # To make the membership vector, you need to provide a named numeric vector,
	# in which the numbers are categories and the names are miRNA or mRNA names.
    
    # make the communities with the graph and membership as inputs.
    a_communities_of_star_graphs.hungarian_algorithm <- make_clusters(graph = secondary_graph,
                                                                      membership = membership_vector)
    
    return(a_communities_of_star_graphs.hungarian_algorithm)
}

