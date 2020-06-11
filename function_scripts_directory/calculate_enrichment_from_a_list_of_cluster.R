library("DOSE")
library("clusterProfiler")
library("org.Hs.eg.db")
library("igraph")
library("foreach")
library("doParallel")
library("GOSemSim")

calculate_enrichment_from_a_cluster_list <- function(the_cluster_list) {

    # parameters to change on different system and machines:
    # extra_libraries_path_directory <- "/net/home/lding/R/x86_64-redhat-linux-gnu-library/3.4/"
    
    # my packages are also installed on the following directories.
    # .libPaths(c(.libPaths(), extra_libraries_path_directory))
    
    #setup parallel backend to use many processors, use all the available cores
    # available_cores_number <- detectCores()
    available_cores_number <- 3 # we only need to iterate three kinds of GO terms.
    current_cluster <- makeCluster(available_cores_number)
    registerDoParallel(current_cluster)
    
    # load packages on each worker, so that all the threads can use packages by parLapply()
    # http://stackoverflow.com/questions/17196261/understanding-the-differences-between-mclapply-and-parlapply-in-r
    clusterEvalQ(current_cluster, library("clusterProfiler"))
    clusterEvalQ(current_cluster, library("org.Hs.eg.db"))
    clusterEvalQ(current_cluster, library("DOSE"))
    clusterEvalQ(current_cluster, library("igraph"))
    
    # send the variable to the environment on the nodes.
    # clusterExport(cl = current_cluster, deparse(substitute(raw3)))
    
    # eventually, you need to stop the threads.
    #stopCluster(current_cluster)
    # -----------------------------------------------------------------------------
    

    # gene-disease associations
    ck_trial <- try(compareCluster(geneCluster = the_cluster_list,
                                   fun = "enrichDGN",
                                   pvalueCutoff  = 0.01,
                                   pAdjustMethod = "BH",
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   qvalueCutoff = 0.05))
    
    if (is.null(attr(ck_trial, "condition")$message)) {
        ck <- ck_trial
        print("ck is done, because there is no error")
    } else {
        ck <- NULL
        print( "ck is not done, because no enrichment of enrichDGN")
    }
    
    
    # Disease-Ontology
    do_trial <- try(compareCluster(geneClusters = the_cluster_list,
                                   fun = "enrichDO",
                                   ont = "DO",
                                   pvalueCutoff  = 0.01,
                                   pAdjustMethod = "BH",
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   qvalueCutoff = 0.05,
                                   readable = FALSE))
    
    if (is.null(attr(do_trial, "condition")$message)) {
        do <- do_trial
        print("do is done, because there is no error")
    } else {
        do <- NULL
        print( "do is not done, because no enrichment of enrichDO")
    }
    
    
    # Network-of-Cancer-Gene
    ncg_trial <- try(compareCluster(geneClusters = the_cluster_list,
                                    fun = "enrichNCG",
                                    pvalueCutoff  = 0.01,
                                    pAdjustMethod = "BH",
                                    minGSSize = 10,
                                    maxGSSize = 500,
                                    qvalueCutoff = 0.05))
    
    if (is.null(attr(ncg_trial, "condition")$message)) {
        # no error message pops up.
        ncg <- ncg_trial
        print("ncg is done, because there is no error.")
    } else {
        ncg <- NULL
        print( "ncg is not done, because no enrichment of enrichDO")
    }
    
    # kegg-over-representation
    kk_trial <- try(compareCluster(geneClusters = the_cluster_list,
                                   fun = "enrichKEGG",
                                   pvalueCutoff  = 0.01,
                                   pAdjustMethod = "BH",
                                   minGSSize = 10,
                                   maxGSSize = 500,
                                   qvalueCutoff = 0.05))
    
    if (is.null(attr(kk_trial, "condition")$message)) {
        # no error message pops up.
        kk <- kk_trial
        print("kk is done, because there is no error.")
    } else {
        kk <- NULL
        print( "kk is not done, because no enrichment of enrichKEGG")
    }

    
    # GO-over-representation
    ontology_words <- c("BP", "CC", "MF")
    # enriched_GO_list <- list(mode="list", length = length(ontology_words))
    
    enriched_GO_list <- foreach(j = seq_along(ontology_words)) %dopar% {
        # for (j in seq_along(ontology_words)) {
        ego_trial <- try(compareCluster(geneClusters = the_cluster_list,
                                        fun = "enrichGO",
                                        OrgDb = org.Hs.eg.db,
                                        ont = ontology_words[j],
                                        pvalueCutoff  = 0.01,
                                        pAdjustMethod = "BH",
                                        minGSSize = 10,
                                        maxGSSize = 500,
                                        qvalueCutoff = 0.05,
                                        readable= FALSE))
        
        if (is.null(attr(ego_trial, "condition")$message)) {
            # no error message pops up.
            ego <- ego_trial
            print("ego is done, because there is no error.")
        } else {
            ego <- NULL
            print( "ego is not done, because no enrichment of enrichGO")
        }
        # enriched_GO_list[[j]] <- ego
        ego
    }
    
    names(enriched_GO_list) <- ontology_words
    ego_BP <- enriched_GO_list[["BP"]]
    ego_CC <- enriched_GO_list[["CC"]]
    ego_MF <- enriched_GO_list[["MF"]]
    
    print("GO is done")
    
    # ------------------compile the biological terms into a list.------
    current_indicator_list <- list(ck, do, ncg, kk, ego_BP, ego_CC, ego_MF)
    names(current_indicator_list) <- c("ck", "do", "ncg", "kk", "ego_BP", "ego_CC", "ego_MF")
    # }
    
    stopCluster(current_cluster)
    
    return(current_indicator_list)
}