#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop("Usage: Rscript make_trait_groups.R <pheno.txt> <geno_numeric.txt> <num_traits> <n_per_group> <output_dir>")
}

pheno_file   <- args[1]
geno_file    <- args[2]
num_traits   <- as.integer(args[3])
n_per_group  <- as.integer(args[4])
output_dir   <- args[5]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("Reading phenotype file...")
pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

message("Reading genotype numeric file...")
geno <- read.table(geno_file, header = TRUE, stringsAsFactors = FALSE)

# Basic sanity checks
if (ncol(pheno) < (num_traits + 1)) {
  stop("Phenotype file has fewer traits than specified.")
}

sample_ids <- pheno[[1]]
geno_samples <- colnames(geno)[4:ncol(geno)]

if (!all(sample_ids %in% geno_samples)) {
  stop("Mismatch between phenotype sample IDs and genotype sample IDs.")
}

for (t in seq_len(num_traits)) {

  trait_name <- paste0("Trait", t)
  message("Processing ", trait_name)

  trait_vals <- pheno[[t + 1]]
  names(trait_vals) <- sample_ids

  # Remove NA samples
  keep <- !is.na(trait_vals)
  trait_vals <- trait_vals[keep]

  if (length(trait_vals) < (3 * n_per_group)) {
    stop(paste("Not enough samples for", trait_name))
  }

  ord <- order(trait_vals)

  low_ids  <- names(trait_vals)[ord][1:n_per_group]
  high_ids <- names(trait_vals)[ord][(length(ord) - n_per_group + 1):length(ord)]

  remaining <- setdiff(names(trait_vals), c(low_ids, high_ids))
  set.seed(123)
  random_ids <- sample(remaining, n_per_group)

  # Extract genotype matrices
  high_mat   <- as.matrix(geno[, high_ids, drop = FALSE])
  low_mat    <- as.matrix(geno[, low_ids, drop = FALSE])
  random_mat <- as.matrix(geno[, random_ids, drop = FALSE])

  # Count ref/alt alleles
  count_alleles <- function(mat) {
    ref <- rowSums(mat == 0, na.rm = TRUE) * 2 + rowSums(mat == 1, na.rm = TRUE)
    alt <- rowSums(mat == 2, na.rm = TRUE) * 2 + rowSums(mat == 1, na.rm = TRUE)
    list(ref = ref, alt = alt)
  }

  high_ct   <- count_alleles(high_mat)
  low_ct    <- count_alleles(low_mat)
  random_ct <- count_alleles(random_mat)

  out <- data.frame(
    snpid = geno$snpid,
    chr   = geno$chr,
    pos   = geno$pos,
    high_ref   = high_ct$ref,
    high_alt   = high_ct$alt,
    low_ref    = low_ct$ref,
    low_alt    = low_ct$alt,
    random_ref = random_ct$ref,
    random_alt = random_ct$alt,
    stringsAsFactors = FALSE
  )

  outfile <- file.path(output_dir, paste0(trait_name, ".txt"))
  write.table(out, outfile, quote = FALSE, sep = "\t", row.names = FALSE)

  message("Written: ", outfile)
}

message("DONE.")
