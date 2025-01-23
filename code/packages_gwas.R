## Install packages and load library for analysis of GWAS data

#install.packages("data.table")    # Ensures version >= 1.9.4
#install.packages("foreach")       # Ensures version >= 1.4.2
#install.packages("stringi")       # Ensures version >= 0.4-1
#install.packages("matrixStats")   # Ensures version >= 0.14.2
#install.packages("Rcpp")  
#install.packages("data_files/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")
library(geneSLOPE)
library(ggplot2)
library(dplyr)
library(devtools)
library(poolSeq) #for CMH test
library(tidyr)
library(stringr)
#install.packages("tidyverse")
library(tidyverse)

### install packages and load library for analysis of TWAS data

install.packages("BiocManager")
BiocManager::install("DESeq2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")
library(DESeq2)
library(grid)
library(ggplot2)
library(pheatmap)
library(sva)
library(tibble)
library(dplyr)
library(tidyr)
