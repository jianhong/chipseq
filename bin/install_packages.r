#!/usr/bin/env Rscript
pkgs <- commandArgs(trailingOnly = TRUE)
if(length(pkgs)>0){
  while(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", 
                     repos = "https://cloud.r-project.org/", 
                     quiet = TRUE)
  }
  pkgs <- pkgs[pkgs %in% BiocManager::available() | grepl("\\/", pkgs)]

  if(any(grepl("\\/", pkgs))){
    pkgs <- c("remotes", pkgs)
  }
  getPkg <- function(pkgs){
    for(pkg in pkgs){
      if(!requireNamespace(pkg, quietly = TRUE)){
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } 
    }
  }
  getPkg(pkgs)
}