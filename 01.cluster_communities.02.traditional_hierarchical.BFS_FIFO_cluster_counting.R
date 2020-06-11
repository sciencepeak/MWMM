rm(list = ls())
ptm <- proc.time()
options(stringsAsFactors = F)


# number of steps the algorithm loop runs is based on the weight threshold.
weight_threshold = "desired"

# traditional hierarchical clustering method applies to integrated mean value weight.
# First step adds the most largest edge with two vertex to the graph.
# Second step adds the second largest edge with two vertex to the graph
# Third step adds the third largest edge with two vertex to the graph
# And so on.... Continue until add the edges are added to the graph.

if (weight_threshold == 0.5) steps = 378
if (weight_threshold == 0.6) steps = 85
if (weight_threshold == 0.7) steps = 38

if (weight_threshold == "desired") steps = 85

# The step_number is a numeric vector as x axis tick marks in plotting
step_number <- seq_len(steps)

weights_directory <- file.path(getwd(),
                               "output_data_directory",
                               "plain_text_files",
                               "six_kinds_of_weights_files.BRCA")

input_file_name_vector <- c("all_negative_value_weight.csv",
                            "all_positive_value_weight.csv",
                            "arithmetic_mean_value_weight.csv",
                            "geometric_mean_value_weight.csv",
                            "integrated_mean_value_weight.csv",
                            "maximum_absolute_value_weight.csv")

for ( current_file_name in input_file_name_vector) {

    # the base name doesn't have the .csv suffix
    current_base_name <- gsub("\\.csv", "", current_file_name)
    current_path_to_file <- file.path(weights_directory, current_file_name)
    
    # input dataframe is a csv edgelist file with columns: microRNA, mRNA, weight, role
    input_df <- read.csv(current_path_to_file, check.names=FALSE)
    
    # sort the edgelist dataframe descreasingly by weight:  (highest -> lowest)
    sorted_edgelist_df <- input_df[(order(input_df$weight, decreasing = T)), ]
    
    backup_sorted_edgelist_df<-sorted_edgelist_df
    
    mRNA_vector <- sorted_edgelist_df$mRNA
    microRNA_vector <- sorted_edgelist_df$microRNA

    # for the current input csv file, if we run from 1 to n step in total,
    # the cluster_num_of_each_step_in_cur_file stores the cluster number when run
    # 1 step, 2 step, 3 step,..., n step, respectively.
    cluster_num_of_each_step_in_cur_file <- vector(mode = "numeric")
    
    # frontier stores the elements to be iterated in the current cluster in the current step in the current file
    # explored stores the elements already iterated in the current cluster in the current step in the current file
    # all_element stores all the mRNA and miRNA of the current step in the current file.
    frontier <- vector(mode = "character")
    explored <- vector(mode = "character") 
    all_element <- vector(mode = "character")
    
    # iterate all the sorted edge list.
    for(current_step in 1:steps){
        
        # define a counter.
        cluster_num_of_cur_step_in_cur_file <- 0
        
        # all_element is like c("hsa-mir-944","hsa-mir-452", "hsa-mir-452", "MRO", "TMEM143", "MYOM3")
        all_element<- append(sorted_edgelist_df[1:current_step, "microRNA"],
                             sorted_edgelist_df[1:current_step, "mRNA"])
        
        # When all the elements have been iterated in the explored, all_element is empty
        # Then we need exit the loop that counts the clusters in the current step 
        while(length(all_element) != 0){
            
            # The first element inserted into frontier is miRNA
            frontier <- append(frontier, all_element[1])
            
            # Keep traverseng in one cluster when the frontier is not empty.
            # In each run of the while loop we hang one depth level of nodes
            # (miRNA or mRNA) to its parental nodes(miRNa or mRNA)
            while(length(frontier)!=0){

                # If the first element in the frontier is miRNA,
                # iteratively expand the miRNA nodes by generating its child leaf nodes of mRNAs
                # linked to that parental node miRNA by edge,
                # and store the generated mRNA into the frontier. 
                if(frontier[1] %in% microRNA_vector){
                    
                    for(i in 1:current_step){
                        if(frontier[1]==sorted_edgelist_df[i, "microRNA"]){
                            if(sorted_edgelist_df[i, "mRNA"] != "HAVE BEEN READ"){
                                frontier <- append(frontier, sorted_edgelist_df[i, "mRNA"])
                            }
                        }
                    }
                    
                    # When all the mRNA nodes of the same depth level are added to the frontier,
                    # rename their parental miRNA node to "HAVE BEEN READ"
                    for(i in 1:current_step){
                        if(frontier[1]==sorted_edgelist_df[i, "microRNA"]){
                            sorted_edgelist_df[i,"microRNA"]="HAVE BEEN READ"
                        }
                    }
                }
                
                # If the first element in the frontier is mRNA,
                # iteratively expand the mRNA nodes by generating its child leaf nodes of miRNAs
                # linked to that parental node mRNA by edge,
                # and store the generated miRNA into the frontier. 
                if(frontier[1] %in% mRNA_vector){
                    for(i in 1:current_step){
                        if(frontier[1]==sorted_edgelist_df[i,"mRNA"]){
                            if(sorted_edgelist_df[i,"microRNA"] != "HAVE BEEN READ"){
                                frontier <- append(frontier,sorted_edgelist_df[i,"microRNA"])
                            }
                        }
                    }
                    
                    # When all the miRNA nodes of the same depth level are added to the frontier,
                    # rename their parental mRNA node to "HAVE BEEN READ"
                    for(i in 1:current_step){
                        if(frontier[1]==sorted_edgelist_df[i, "mRNA"]){
                            sorted_edgelist_df[i, "mRNA"]="HAVE BEEN READ"
                        }
                    }
                }
                
                # In each run of the current whileloop,
                # the first element of the frontier has been iterated,
                # move the element from frontier(a FIFO queque) to the explored
                explored <- append(explored, frontier[1])
                frontier <- frontier[-1]
            }
            
            # When iterate all the nodes in one cluster,
            # update the cluster count for the current step in the current file.
            cluster_num_of_cur_step_in_cur_file <- cluster_num_of_cur_step_in_cur_file + 1
            
            # remove the explored elements from the all_element
            all_element <- setdiff(all_element, explored)
        }
        
        # When we go through all the nodes from one cluster,
        # initialize the explored for the next cluster
        explored <- vector(mode = "character")
        
        # sorted_edgelist_df also need to be initialized from the backup untouched data frame. 
        # because its elements are changed to "HAVE BEEN READ"
        sorted_edgelist_df <- backup_sorted_edgelist_df
        
        cluster_num_of_each_step_in_cur_file <- append(cluster_num_of_each_step_in_cur_file,
                                                       cluster_num_of_cur_step_in_cur_file)
        
        print(cluster_num_of_cur_step_in_cur_file)
        
    }

    # wrap the step and cluster numbers into data frames for plotting
    if(current_file_name =="all_negative_value_weight.csv"){
        graph1 <- data.frame(step_number,cluster_num_of_each_step_in_cur_file)
    }
    if(current_file_name =="all_positive_value_weight.csv"){
        graph2 <- data.frame(step_number,cluster_num_of_each_step_in_cur_file)
    }
    if(current_file_name =="arithmetic_mean_value_weight.csv"){
        graph3 <- data.frame(step_number,cluster_num_of_each_step_in_cur_file)
    }
    if(current_file_name =="geometric_mean_value_weight.csv"){
        graph4 <- data.frame(step_number,cluster_num_of_each_step_in_cur_file)
    }
    if(current_file_name =="integrated_mean_value_weight.csv"){
        graph5 <- data.frame(step_number,cluster_num_of_each_step_in_cur_file)
    }
    if(current_file_name =="maximum_absolute_value_weight.csv"){
        graph6 <- data.frame(step_number,cluster_num_of_each_step_in_cur_file)
    }
    
    # zero the cluster number vector, we are ready for looping the next file.
    cluster_num_of_each_step_in_cur_file <- vector(mode = "numeric")
}


plot(graph1$step_number,bty="l",graph1$cluster_num_of_each_step_in_cur_file,col="red",xlim=c(1,max(graph1$step_number)),ylim=c(1,20),lwd=2,xlab="Steps",ylab="Numbers of connected clusters",type="b",pch=15,cex=0.4)

lines(graph2$step_number,bty="l",graph2$cluster_num_of_each_step_in_cur_file+0.08,col="black",xlim=c(1,max(graph1$step_number)),ylim=c(1,20),lwd=2,xlab="Steps",ylab="Numbers of connected clusters",type="b",pch=16,cex=0.4)

lines(graph3$step_number,bty="l",graph3$cluster_num_of_each_step_in_cur_file+0.16,col="orange",xlim=c(1,max(graph1$step_number)),ylim=c(1,20),lwd=2,xlab="Steps",ylab="Numbers of connected clusters",type="b",pch=17,cex=0.4)

lines(graph4$step_number,bty="l",graph4$cluster_num_of_each_step_in_cur_file+0.24,col="blue",xlim=c(1,max(graph1$step_number)),ylim=c(1,20),lwd=2,xlab="Steps",ylab="Numbers of connected clusters",type="b",pch=18,cex=0.4)

lines(graph5$step_number,bty="l",graph5$cluster_num_of_each_step_in_cur_file+0.32,col="green",xlim=c(1,max(graph1$step_number)),ylim=c(1,20),lwd=2,xlab="Steps",ylab="Numbers of connected clusters",type="b",pch=19,cex=0.4)

lines(graph6$step_number,bty="l",graph6$cluster_num_of_each_step_in_cur_file+0.40,col="pink",xlim=c(1,max(graph1$step_number)),ylim=c(1,20),lwd=2,xlab="Steps",ylab="Numbers of connected clusters",type="b",pch=24,cex=0.4)


legend("topleft",bty="n",text.width=1.5,cex=0.7,legend=c("all_negative_value_weight  ",
                                                       "all_positive_value_weight  ",
                                                       "arithmetic_mean_value_weight",
                                                       "geometric_mean_value_weight  ",
                                                       "integrated_mean_value_weight ",
                                                       "maximum_absolute_value_weight"),
       
       col=c("red","black","orange","blue","green","pink"),seg.len=0.3,x.intersp=0.5,lwd=4)

proc.time() - ptm
