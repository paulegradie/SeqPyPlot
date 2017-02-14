#! /usr/bin/Rscript

#Read in count matrix and normlize using the edgeR package

setwd(".")

#Set Args
args = commandArgs(trailingOnly=TRUE)



#First block modified from:
# Script name: instant_pkgs.R
# Purpose: Package installation and loading
# Author: Kay Cichini
# Date: 2012-06-19
# Licence: cc by-nc-sa

instant_pkgs <- function(pkgs) { 
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    
    #install.packages(pkgs_miss)
    
    result <- try(source("https://bioconductor.org/biocLite.R"));
    
    if(class(result) == "try-error")
      {source("http://bioconductor.org/biocLite.R")}

    ## try http:// if https:// URLs are not supported
    biocLite(pkgs_miss)

    }
  
  if (length(pkgs_miss) == 0) {
    message("\n ...R Packages were already installed!\n")
  }
  
  # install packages not already loaded:
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss)
  }
  
  # load packages not already loaded:
  attached <- search()
  attached_pkgs <- attached[grepl("package", attached)]
  need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
  
  if (length(need_to_attach) > 0) {
    for (i in 1:length(need_to_attach)) require(need_to_attach[i], character.only = TRUE)
  }
  
  if (length(need_to_attach) == 0) {
    message("\n ...Packages were already loaded!\n")
  }
}

#Install edgeR - and load it
instant_pkgs(c("edgeR"))
require("edgeR")

#####################
#PERFORM NORMALIZATION
data_matrix=as.data.frame(read.delim(args[1], sep="\t", row.names="Gene"))
#data_matrix=as.matrix(read.delim("current_test_count_matrix.txt", sep="\t", row.names="Gene"))

#dat_matrix=as.matrix(dat)
to_normalize_matrix = DGEList(data_matrix, remove.zeros=TRUE)
normalized_matrix = calcNormFactors(to_normalize_matrix)

write.table(as.matrix(normalized_matrix), file="normalized_count_data.txt", sep = '\t', col.names=FALSE)
