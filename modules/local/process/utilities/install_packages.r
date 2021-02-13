#!/usr/bin/env Rscript
pkgs <- commandArgs(trailingOnly = TRUE)
# set libPath to pwd
lib <- .libPaths()
if(file.access(lib[1], mode=2)!=0){
  pwd <- getwd()
  pwd <- file.path(pwd, "lib")
  dir.create(pwd)
  .libPaths(c(pwd, lib))
}

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