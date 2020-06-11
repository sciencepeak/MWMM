library("igraph")

merge_blossom_matching_clusters <- function(a_communities, blossom_matching_result_file, raw_edge_list_df) {
    
    blossom_matching_df <- read.csv(blossom_matching_result_file, check.names = F)

    all_node_members <- numeric(0)
    
    for (i in seq_len(nrow(blossom_matching_df))) {
        cluster_from <- communities(a_communities)[[blossom_matching_df[i, "from"]]]
        cluster_to <- communities(a_communities)[[blossom_matching_df[i, "to"]]]
        merged_cluster_nodes <- c(cluster_from, cluster_to)
        
        # Due to OCD, the miRNA names are put to the end of the cluster.
        cluster_mRNAs <- merged_cluster_nodes[!grepl("hsa-", merged_cluster_nodes)]
        cluster_miRNAs <- merged_cluster_nodes[grepl("hsa-", merged_cluster_nodes)]
        sorted_merged_cluster_nodes <- c(cluster_mRNAs, cluster_miRNAs)
        
        # The new way of organizing clusters is to construct a communities object.
        # The hardest job is to create a named vector as the membership of a communities object.
        current_numbered_members <- rep(i, times = length(sorted_merged_cluster_nodes))
        names(current_numbered_members) <- sorted_merged_cluster_nodes
        
        all_node_members <- append(all_node_members, current_numbered_members)
    }
    
    # From the raw graph subset the subgraph to construct communities.
    raw_graph = graph.data.frame(raw_edge_list_df, directed = F)
    raw_adjacency_matrix = as_adjacency_matrix(raw_graph, type="both", names=TRUE, sparse=FALSE, attr="weight")
    
    sub_adjacency_matrix = raw_adjacency_matrix[names(all_node_members), names(all_node_members)]
    sub_graph = graph_from_adjacency_matrix(sub_adjacency_matrix, mode = "undirected", weighted = T)
    
    a_communities_of_merged_clusters <- make_clusters(graph = sub_graph, membership = all_node_members)
    
    return(a_communities_of_merged_clusters)
}
