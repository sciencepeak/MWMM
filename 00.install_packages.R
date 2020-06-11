install.packages("igraph") # Louvain and other clustering algorithms are provided.
install.packages("hash")
install.packages("magrittr")
install.packages("foreach")
install.packages("doParallel")
install.packages("tidyr")
install.packages("reshape2")

install.packages("Matrix")
install.packages("MASS")
install.packages("tweenr")
install.packages("units")
install.packages("lattice")
install.packages("formattable") # This package provides percent function.


install.packages("ggplot2")
install.packages("GGally") # ggnet2 is available through the GGally package.
install.packages("network") # ggnet2 depends on this package.
install.packages("sna")  # ggnet2 depends on this package.
install.packages("clue") # for hungarian algorithm.


# With R version 3.5 or greater, install Bioconductor packages using BiocManager; see https://bioconductor.org/install
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install("clusterProfiler")
# etc...

source("https://bioconductor.org/biocLite.R")
biocLite("GOSemSim")
biocLite("org.Hs.eg.db")
biocLite("DOSE")
biocLite("clusterProfiler")
biocLite("rvcheck")

# library("rvcheck")
# rvcheck::check_bioc("GOSemSim")
# rvcheck::check_bioc("org.Hs.eg.db")