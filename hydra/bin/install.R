install.packages("Rcpp", repos='http://cran.us.r-project.org')
install.packages('tidyverse', repos='http://cran.us.r-project.org')
install.packages('dplyr', repos='http://cran.us.r-project.org')

# Required for survminer
install.packages("ggplot2", repos='http://cran.us.r-project.org')
devtools::install_url("https://github.com/wilkelab/cowplot/archive/0.9.0.zip")

# Install xCell
devtools::install_github('dviraran/xCell')
install.packages("survminer", repos='http://cran.us.r-project.org')
install.packages('IRkernel', repos='http://cran.us.r-project.org')
install.packages('cluster', repos='http://cran.us.r-project.org')
install.packages('factoextra', repos='http://cran.us.r-project.org')
#install.packages('NB', repos='http://cran.us.r-project.org')

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("annotate")
BiocManager::install("geneplotter")
BiocManager::install("GSEABase")
BiocManager::install("ConsensusClusterPlus")
BiocManager::install("M3C")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("eisa")
BiocManager::install("fgsea")
BiocManager::install("clusterProfiler")
