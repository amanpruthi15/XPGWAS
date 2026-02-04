#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    "Usage: Rscript run_xpgwas.R <xpgwas_input_dir> <trait_name> <p_threshold> <output_dir>",
    call. = FALSE
  )
}

input_dir  <- args[1]
trait_name <- args[2]
p_thresh   <- as.numeric(args[3])
output_dir <- args[4]

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

library(XPGWAS)

input_files <- list.files(
  input_dir,
  pattern = paste0(trait_name, ".*\\.rds$"),
  full.names = TRUE
)

if (length(input_files) == 0) {
  stop("No XP-GWAS input files found for the specified trait")
}

results_list <- list()

for (f in input_files) {
  message("Processing: ", basename(f))
  chunk_data <- readRDS(f)

  res <- xpgwas(
    geno = chunk_data$geno,
    pheno = chunk_data$pheno,
    trait = trait_name
  )

  results_list[[basename(f)]] <- res
}

saveRDS(
  results_list,
  file = file.path(
    output_dir,
    paste0("XP_GWAS_results_", trait_name, ".rds")
  )
)

message("XP-GWAS completed successfully.")
