library("GGally") # ggnet2 is available through the GGally package.
library("network") # ggnet2 depends on this package.
library("sna")  # ggnet2 depends on this package.
library("ggplot2") # ggnet2 depends on this package.
library("igraph")
library("reshape2") # provide melt function.

# The most important input of the graph-drawing function is the edge lists.
# path_to_graph_output includes the output file name with file type extension.
draw_cluster <- function(cluster_edge_list_df, raw_edge_list_df, path_to_graph_output, file_type = "pdf") {
    
    cluster_graph <- graph.data.frame(cluster_edge_list_df, directed = F)
    
    cluster_graph_adjacency_matrix <- as_adjacency_matrix(cluster_graph, type="both", names=TRUE, sparse=FALSE, attr="weight")
    
    raw_green_edge_list = subset(raw_edge_list_df, edge_color == 'green')
    
    # edge.label must get from mirror edgelist
    mirror_origin <- subset(melt(cluster_graph_adjacency_matrix), value > 0)
    
    # set pramaters for graph
    color_node = ifelse(rownames(cluster_graph_adjacency_matrix) %in% unique(raw_edge_list_df$microRNA), "steelblue", "orange")
    shape_node = ifelse(rownames(cluster_graph_adjacency_matrix) %in% unique(raw_edge_list_df$microRNA), 15, 16)
    color_edge = ifelse(mirror_origin$value %in% raw_green_edge_list$weight, "green","red")
    label_edge =  round(mirror_origin$value, 2)
    
    
    # plot the bipartite network object in smaller font
    ggnet2(cluster_graph_adjacency_matrix,
           label = T,
           label.size = 2.5,
           # legend.position = "null",
           # size = "mode",
           # size.palette = c("actor" = 3, "event" = 5),
           max_size = 6,
           color = color_node,
           shape = shape_node,
           mode = "fruchtermanreingold",
           edge.color = color_edge,
           edge.label = label_edge,
           edge.size = 0.5,
           edge.label.size = 2.5
    )
    ggsave(path_to_graph_output, plot = last_plot(),
           width = 170, height = 170, units = "mm", device = file_type)
    
}