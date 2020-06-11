library("igraph")
library("magrittr") # pipeline operation.

calculate_inner_weight_from_a_communities <- function(a_communities, reference_edge_list_df) {
    
    get_inner_weight_from_a_cluster_nodes <- function(a_cluster_nodes){
        
        # a_cluster_nodes is a character vector of mRNA and mRNA nodes
        # a_cluster_nodes also represent a community/cluster in the communities object.
        a_cluster_miRNAs <- a_cluster_nodes[grepl("hsa-", a_cluster_nodes)]
        a_cluster_mRNAs <- a_cluster_nodes[!grepl("hsa-", a_cluster_nodes)]
        
        condition_miRNA <- reference_edge_list_df$microRNA %in% a_cluster_miRNAs
        condition_mRNA <- reference_edge_list_df$mRNA %in% a_cluster_mRNAs
        condition_both <- condition_miRNA & condition_mRNA
        
        a_cluster_edge_list_df <- reference_edge_list_df[which(condition_both), ]
        
        a_cluster_inner_weight <- a_cluster_edge_list_df %>%
            .[, "weight"] %>%
            sum
    }
    
    # inner_weight (IW) emitter to the receivers within one cluster.
    IW_vector <- sapply(communities(a_communities), get_inner_weight_from_a_cluster_nodes)
    return(IW_vector)
}


calculate_outer_weight_from_a_communities <- function(a_communities, reference_edge_list_df) {
    
    E2MROW_vector <- numeric(0)
    R2MEOW_vector <- numeric(0)
    
    for ( i in seq_along(communities(a_communities))) {
        
        # Get the vertexes (nodes) of the current cluster and other clusters.
        current_cluster_nodes <- communities(a_communities)[[i]]
        all_other_cluster_communities <- communities(a_communities)[-i]
        all_other_cluster_nodes <- unname(unlist(all_other_cluster_communities))
        
        # Differentiate miRNA and mRNA vertexes by string "hsa-"
        current_cluster_miRNAs <- current_cluster_nodes[grepl("hsa-", current_cluster_nodes)]
        current_cluster_mRNAs <- current_cluster_nodes[!grepl("hsa-", current_cluster_nodes)]
        all_other_cluster_miRNAs <- all_other_cluster_nodes[grepl("hsa-", all_other_cluster_nodes)]
        all_other_cluster_mRNAs <- all_other_cluster_nodes[!grepl("hsa-", all_other_cluster_nodes)]
        
        # Now, calculate two types of outer weight.
        
        # emitter to matched receiver outer weight: E2MROW
        # current_cluster_miRNAs vs. all_other_cluster_mRNAs
        
        condition_miRNA <- reference_edge_list_df$microRNA %in% current_cluster_miRNAs
        condition_mRNA <- reference_edge_list_df$mRNA %in% all_other_cluster_mRNAs
        condition_both <- condition_miRNA & condition_mRNA
        
        current_cluster_E2MROW_weight <- reference_edge_list_df %>%
            .[which(condition_both), ] %>%
            .[, "weight"] %>%
            sum
        
        E2MROW_vector <- append(E2MROW_vector, current_cluster_E2MROW_weight)
        
        # receiver to matched emitter outer weight: R2MEOW
        # current_cluster_mRNAs vs. all_other_cluster_miRNAs
        
        condition_miRNA <- reference_edge_list_df$microRNA %in% all_other_cluster_miRNAs
        condition_mRNA <- reference_edge_list_df$mRNA %in% current_cluster_mRNAs
        condition_both <- condition_miRNA & condition_mRNA
        
        current_cluster_R2MEOW_weight <- reference_edge_list_df %>%
            .[which(condition_both), ] %>%
            .[, "weight"] %>%
            sum
        
        R2MEOW_vector <- append(R2MEOW_vector, current_cluster_R2MEOW_weight)
    }
    # -------------------------------------------------------------------
    
    outer_weights_df <- data.frame(E2MROW = E2MROW_vector,
                                   R2MEOW = R2MEOW_vector)
    return(outer_weights_df)
}

# piece together the inner weight and outer weight and calculate condition 1 and 2.
get_inner_outer_weights_and_conditions <- function(a_communities, reference_edge_list_df) {
    
    # get three vectors of different weights.
    IW_vector <- calculate_inner_weight_from_a_communities(a_communities, reference_edge_list_df)
    outer_weight_dataframe <- calculate_outer_weight_from_a_communities(a_communities, reference_edge_list_df)
    E2MROW_vector <- outer_weight_dataframe$E2MROW
    R2MEOW_vector <- outer_weight_dataframe$R2MEOW
    
    # Condition one specifies IW > E2MROW.
    # Condition two specifies 2 × IW > E2MROW + R2MEOW
    condition_01 <- IW_vector > E2MROW_vector
    condition_02 <- 2 * IW_vector > E2MROW_vector + R2MEOW_vector
    
    # wrap up the three weights and two conditions into a dataframe
    three_weights_dataframe <- data.frame(IW = IW_vector,
                                          E2MROW = E2MROW_vector,
                                          R2MEOW = R2MEOW_vector,
                                          condition_01 = condition_01,
                                          condition_02 = condition_02)
    return(three_weights_dataframe)
}
