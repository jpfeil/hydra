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
install.packages('NB', repos='http://cran.us.r-project.org')

source("https://bioconductor.org/biocLite.R")
biocLite("annotate")
biocLite("geneplotter")
biocLite("GSEABase")
biocLite("ConsensusClusterPlus")
biocLite("M3C")
biocLite("org.Hs.eg.db")
biocLite("DOSE")
biocLite("eisa")
biocLite("fgsea")
biocLite("clusterProfiler")
