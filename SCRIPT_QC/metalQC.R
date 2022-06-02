#!/usr/bin/Rscript

# Loading data from stding
args = commandArgs(trailingOnly=TRUE)
metal = fread(args[1])
N = as.numeric(args[2])
## Library loading
.libPaths("/mnt/rreal/RDS/acarrasco/R_libs/")
library(data.table)
library(stringr)
library(tidyverse)

# Loading all the functions needed to do the QC
# Using here as these functions are under the root project folder
#source(here::here("R/Utils.R"))
source("/mnt/rreal/RDS/acarrasco/TOOLS/postgwasQC/R/Utils.R")


harmonise_metal(metal, N)
