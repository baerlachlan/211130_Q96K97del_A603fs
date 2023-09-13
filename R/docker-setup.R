## This script was made to reproduce results that were impacted by a package version upgrade
## It is intended to be run on a fresh docker container spun up from the bioconductor/bioconductor_docker:RELEASE_3_16 image
## After successful installation of packages, the R/enrichment.Rmd file can be rendered to produce the results presented in the associated manuscript

library(devtools)

BiocManager::install(c(
  "tidyverse", "magrittr", "future.apply", "here", "AnnotationHub", "purr",
  "scales", "kableExtra", "tictoc", "ggrepel", "RColorBrewer", "ggpubr",
  "pander", "rmarkdown", "edgeR", "limma", "pheatmap", "cqn", "DT", "htmltools",
  "plyranges", "ensembldb", "fgsea", "statmod"
))
install_version("msigdbr", version = "7.4.1")
install_version("babelgene", version = "22.3")
