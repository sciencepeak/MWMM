rm(list = ls())
options(stringsAsFactors = F)
library("igraph")

# take an example of edge list.
input_edge_list <- data.frame(X = c("a", "a", "b", "c", "d", "d"),
                              Y = c("B", "C", "C", "A", "C", "A"),
                              weight = c(1, 2, 3, 4, 5, 7))
input_edge_list

# From edge list to adjacency matrix.
a_graph <- graph.data.frame(input_edge_list, directed = F)
plot(a_graph)
an_adjacency_matrix = as_adjacency_matrix(a_graph, type="both", names=TRUE, sparse=FALSE, attr="weight")
an_adjacency_matrix

# From adjacency matrix to edge list:
the_graph <- graph_from_adjacency_matrix(an_adjacency_matrix, mode = "undirected", weighted = T)
plot(the_graph)
new_edge_list <- as_data_frame(the_graph)
names(new_edge_list) <- c("X", "Y", "weight")
new_edge_list

