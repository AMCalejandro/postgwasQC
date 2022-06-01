#!/usr/bin/Rscript 

## Library loading
.libPaths("/mnt/rreal/RDS/acarrasco/R_libs/")
library(data.table)
library(stringr)
library(tidyverse)
library(colochelpR)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(here)
library(purrr)
library(BSgenome.Hsapiens.NCBI.GRCh38) # Necessary if I want to deduce the genome build


# TODO
# Add a proper example usage on the script.
# The user is expected to pass:
# Mandatory -> gwasFile, Nsamples, OUTCOMEcolname
# Optional mafFile
# if length(args == 0) I run an example usage function

args = commandArgs(trailingOnly=TRUE)
lengthArgs = length(args)
if (length(args) < 4) {
  message("\nCheck all mandatory arguments are present: 1.GWAS 2.Nsamples 3. name OUTCOME col \n")
  message("If MAF file is missing, pass NULL as the argument")
  message("\n**Example usage...**\n")
  message("Rscript harmonization.R [path_GWAS] [path_maf or NULL if mising] [Nsamples] [OUTCOME colname]***\n")
  stop("Exiting....\n")
}


# Loading data from stding
gwas = fread(args[1])
maf = fread(args[2]) # This argument is optional
Nsamples = args[3]
outcome = args[4]

# Loading all the functions needed to do the QC
# Using here as these functions are under the root project folder
source(here::here("R/Utils.R"))


gwas_QC = harmonise_gwas(gwas)
maf_QC = harmonise_maf(maf)
res = harmonization(gwas = gwas_QC, 
                    maf = maf_QC,
                    N = Nsamples,
                    multiple_outcome = T,
                    outcome_var = outcome,
                    writeOut = T)