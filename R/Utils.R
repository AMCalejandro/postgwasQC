# TODO
# I need to come up with something to save the directory to save results


# TODO
# In {harmonise_gwas() }, when I perform the conversion from loc to rs
# There is a potential disruption in harmonization
# If input data has snp in loc format under hg38 build, when I harmonise
# The resulting CHR BP columns will still be in hg38, whereas the idea
# is to harmonise all columns under hg19 build
# This is mostly an aesthtetic issue because, we use rsIDs dowstream when we meta analyse
# Still want to improve it



# TODO
# I need to add util to remove rows with NAs in the SNP metrics columns
# # I want to perform an internal check as well.
# Basically compare the initial sample size with the final sample size.
# If we detect a high decrease in the resulting number of SNPs, raise a warning




# TODO
# Find a method to figure out the genome build from the input data
# Strategy. Use a while loop to get go over each snp on a gwas df.
# For each snp, we use   library(BSgenome.Hsapiens.NCBI.GRCh38) data = snpsByOverlaps(dbSNP, "1:752895-752895", genome = "GRCh38")
# If the resulting data is empty, we keep iterating until we populate it.
# If is not null, we then check that the reference and alternate alleles match
# Ideally, we would do this x times to cross validate the genome build inferance


#' Function to determine the genome build from input data
#'
#' @param df df. Input GWAS with SNP column already splitted with harmonise_GWAS
#' @param psychencode_SNP_info df. SNP information for all QTLs considered,
#'   including rsIDs (if available), location, and reference and alternate
#'   alleles, as downloaded from \url{http://resource.psychencode.org/}.
#' @param add_colnames logical. Should columne names be added to the file?
#'   Default is TRUE i.e. assumes user has not derived column names for full
#'   summary statistic from another psychencode-derived file.
#'
#' @return genome build detected as chr
#' @export
#'


# I found indels to be proble
inferBuild = function(df) {
  n = 0
  index = 0
  count = 0
  dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38
  
  while (n < 60) {
    index = index + 1
  
    myRow = df[index, ]
    pattern = paste0(myRow$CHR, ":", myRow$BP, "-", myRow$BP)
    refsnp = snpsByOverlaps(dbSNP,  pattern , genome = "GRCh38") %>% 
      as.data.frame
    
    if (nrow(refsnp) == 0) {
      next
    }
    else {
      n = n + 1
      if ( any((myRow$REF == refsnp$ref_allele) & 
           (myRow$ALT %in% refsnp$alt_alleles)) ) {
        count = count + 1
      }
      else {
        count = count - 1
      }
    }
  }
  if (count == -60) {
    build = "hg19"
  } else if (count == 60) {
    build = "hg38"
  } else {
    
    warning("Some discrepancies between the REF and ALT alleles from GWAS files have been found compared to the HSapiens Genome Sequence by NCBI \n")
    
    if (count > 10) {
      build = "hg38"
    } else if (count < -10) {
      build = "hg19"
    } else {
      stop("Something is going wrong")
    }
    message("Using the most likely build: ", build)
  }
  return(build)
}



# TODO
# In this file, function for post metal qc?

# Using the build info to harmonose using colochelpR
## TODO
# Fin a method that use any function such as find overlap to infet the build
# of the input data
harmonise_metal = function(df, N) {
  data_filtered <- df 
    filter(HetDf >= N -2) 
    
  data_filtered_sorted <- data_filtered %>%
    arrange(`P-value`)
  
  #Filter out SNPs with HetPVal < 0.05 (Cochran's Q-test for heterogeneity)
  #Also filter out SNPs with HetISq > 80
  data_filtered_sorted_het <-data_filtered_sorted %>%
    filter(HetPVal > 0.05) %>%
    filter(HetISq < 80)
  
  #Check MAF variability - remove variants with MAF variability > 15%
  data_filtered_sorted_het_MAF <- data_filtered_sorted_het %>%
    mutate(MAF_variability = MaxFreq - MinFreq) %>%
    filter(MAF_variability <= 0.15)
  
  #Export for FUMA
  export_FUMA <- data_filtered_sorted_het_MAF %>%
    select(MarkerName, `P-value`, Allele1, Allele2, Effect, StdErr, TotalSampleSize) %>%
    rename(rsID = MarkerName,
           pval = `P-value`)
  
  fwrite(export_FUMA, "metaanalysis_FUMA.txt", quote = F, row.names = F, col.names = T, sep = "\t")
}



# Probably I do not want to create a function just for this. 
# I should consider including this somewhere
harmonise_gwas = function(df) {
  
  # Get SNP metadata
  rndn_snp = df[1,]$SNP
  singleSep = sepCharacters(rndn_snp)
  indels = indelsFinder(df$SNP)
  snpFormat = checkFormat(rndn_snp)
  
  
  # Start SNP processing
  if (!snpFormat == "rs") {
    data_tmp = snpProcessor(df, sep = singleSep, ind = indels, snpFormat = snpFormat)
    
    # Minor processing
    if ( !is.numeric(data_tmp$CHR) | !is.numeric(data_tmp$BP) ) { # Minor processing 
      data_tmp = data_tmp %>%
        mutate(CHR = as.numeric(gsub(".*chr", "", CHR)),
               BP = as.numeric(BP)) 
    }
    
    # Infer chromosome build
    build = inferBuild(data_tmp)
    
    # Get rsIDs
    if (build == "hg19") {
      dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37
    } else {
      dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38
    }
    df_final = colochelpR::convert_loc_to_rs(df = data_tmp, dbSNP = dbSNP) %>%
      dplyr::filter(!is.na(SNP)) %>%
      dplyr::relocate(SNP, .before = CHR) %>%
      dplyr::mutate(CHR =  as.numeric(levels(CHR)[CHR]))
    # TODO
    # If the build is hg38, I want to get CHR and BP in hg19
  }
  
  else if (snpFormat == "rs") {
    dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37
    df_final =  colochelpR::convert_rs_to_loc(df = df, SNP_column = "SNP", dbSNP = dbSNP)
  }
  else {
    stop("Make sure the SNP column is named in the correct way and SNP data is in the correct format")
  }
  
  return(df_final)
}




# Function to get MAF information
harmonise_maf = function(df) {
  
  rndn_snp = df[1,]$SNP  
  # PART 1
  # Understand the format the SNP is in
  indels = indelsFinder(df$SNP)
  snpFormat = checkFormat(rndn_snp)
  
  
  if (!snpFormat == "rs") {
    df = df %>% dplyr::select(-CHR)
    
    # Get the special character
    snp_sep = substr(snpFormat, 4,4)
    snp_splitted = 
      as.data.frame(stringr::str_split_fixed(df$SNP, n = 4, pattern= snp_sep))
    names(snp_splitted) = c("CHR", "BP", "REF", "ALT")
    
    data_final = cbind(snp_splitted, df) %>% 
      dplyr::select(-SNP) 
    if ( !is.numeric(data_final$CHR) | !is.numeric(data_final$BP) ) {  
      data_final = data_final %>%
        mutate(CHR = as.numeric(gsub(".*chr", "", CHR)),
               BP = as.numeric(BP))
    }
    
    data_final = data_final %>%
      dplyr::select(-c(REF,ALT))
    
  }
  else if (snpFormat =="rs") {
    data_final = df
    warning("This maf file does not need processing and is ready to be merged by rsID")
  }
  return(data_final)
}





harmonization = function(gwas,
                         maf = NULL,
                         N = Nsamples,
                         multiple_outcome = FALSE,
                         outcome_var = NULL,
                         writeOut = FALSE) {
  
  if (is.null(maf)) {
    data_final = getMAF(gwas)
  } 
  else {
    data_final = gwas  %>%
      dplyr::inner_join(maf, by = c("CHR", "BP", "A1"))
  }
  
  data_final = data_final %>%
    dplyr::relocate(c(A2,MAF), .after = A1) %>%
    dplyr::mutate(Nsamples = N)
  
  if (multiple_outcome) {
    data_final = data_final %>% group_by(!!! rlang::syms(outcome_var))
    data_final = data_final %>% group_split()
    if (writeOut) {
      purrr::map(data_final, ~writeOutput(.x, outcome_var))
    }
  }
  else {
    if (writeOut) {
      writeOutput(data_final)
    }
  }
  return(data_final)
}

# Function to write the output. This function searches for an outcome variable
# on the input datafra,e, and uses the unique value of it to write the name 
# of the final data frame
writeOutput = function(df, outcome_var) {
  getOutcome = dplyr::distinct(df, !!! rlang::syms(outcome_var)) %>% 
    pull()
  
  fwrite(df, paste0("GWAS_", getOutcome, ".txt"),
         quote=F, sep ="\t", col.names =T, row.names = F) 
  
  cat("GWAS for outcome ", getOutcome, " successfully written \n")
}




getMAF = function(df) {
  mafdb <- MafDb.1Kgenomes.phase3.hs37d5
  rsids = df$SNP
  
  mafs <- GenomicScores::gscores(x = mafdb, ranges = rsids %>% as.character(), pop = "EUR_AF")
  
  mafs <- mafs %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "SNP") %>%
    dplyr::rename(maf = EUR_AF) %>% 
    dplyr::select(SNP, maf)
  
  df <- df %>% inner_join(mafs, by = c("SNP"))
  return(df)
}


# Figuring out the snp format from the input document
checkFormat = function(snp) { 
  if (grepl(paste0("chr", 1:22, ":", collapse="|"), x = snp)) {
    snpFormat = "chr:bp"
  } else if (grepl(paste0(1:22, ":", collapse="|"), x = snp)) {
    snpFormat = "chr:bp"
  } else if (grepl(paste0("chr", 1:22, "_", collapse="|"), x = snp)) {
    snpFormat = "chr_bp" 
  } else if (startsWith(snp, prefix= "rs")) {
    snpFormat = "rs"
  } else {
    message("Make sure the variant columnd is named \"SNP\"") 
    message("In addition, only chr:bp, chr_bp, or rs formats are supported")
    stop("Make sure input is in the right format")
  }
  return(snpFormat)
}

# Understand the SNP data a bit further.
# Sometimes the "*" pattern is present, which stands for an indel
# If this character is present, using separate from tidyverse causes disruptions
indelsFinder = function(snpVector) {
  if (any(grepl("*",  snpVector))) {
    indelStatus = TRUE 
  } else {
    indelStatus = FALSE
  }
  return(indelStatus)
}


# This function should be run always that snpFormat = chr:bp
sepCharacters = function(rndn_snp) {
  if ( grepl("_", rndn_snp) & grepl(":", rndn_snp) ) {
    singleSep = FALSE
  } else {
    singleSep = TRUE
  }
  return(singleSep)
}


# This function is meant to be used for data with SNP on CHRBP format
# Unnecessary processing for data with rsIDs
snpProcessor = function(df,
                        sep = TRUE,
                        ind = indels, 
                        snpFormat = snpFormat) {
  if (!snpFormat == "rs") {
    if (!sep) { # In case I have more than one special character for sepparation
      snp_splitted = as.data.frame(stringr::str_split_fixed(df$SNP, pattern = "_",  n = 2))
      snp_splitted = 
        cbind(as.data.frame(stringr::str_split_fixed(snp_splitted$V1, pattern = ":",  n = 4)), 
              snp_splitted$V2)
      names(snp_splitted) = c("CHR", "BP", "REF", "ALT", "A1")
      
    }
    else {
      # Get the special character
      snp_sep = substr(snpFormat, 4,4)
      snp_splitted = 
        as.data.frame(stringr::str_split_fixed(df$SNP, n = 5, pattern= snp_sep))
      names(snp_splitted) = c("CHR", "BP", "REF", "ALT", "A1")

    }
  } 
  else { # If the snp is in rsID format, the sepparator will always be "_"
    snp_sep = "_" 
    snp_splitted = 
      as.data.frame(stringr::str_split_fixed(df$SNP, n = 2, pattern = snp_sep))
    
    names(snp_splitted) = c("SNP", "A1")
    
  }
  df = cbind(snp_splitted, df %>% dplyr::select(-SNP)) 
  return(df)
}
