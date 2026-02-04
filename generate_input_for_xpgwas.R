#!/usr/bin/env Rscript
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript script.R <pheno> <geno> <num_traits> <n_per_group> <outdir> <chunk_size>")
}

pheno_file  <- args[1]
geno_file   <- args[2]
num_traits  <- as.integer(args[3])
n_group     <- as.integer(args[4])
outdir      <- args[5]
chunk_size  <- as.integer(args[6])

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

### Load phenotype
pheno <- fread(pheno_file)
setnames(pheno, 1, "ID")

geno_ids <- pheno$ID

### Pre-read genotype header only
geno_header <- names(fread(geno_file, nrows = 0))
snp_ids <- geno_header[-1]
n_snps <- length(snp_ids)

message("Total SNPs detected: ", n_snps)

### Trait loop
for (t in seq_len(num_traits)) {

  trait <- paste0("Trait", t)
  message("Processing ", trait)

  vals <- pheno[[t + 1]]
  names(vals) <- pheno$ID
  vals <- vals[!is.na(vals)]

  ord <- order(vals)
  low_ids  <- names(vals)[ord][1:n_group]
  high_ids <- names(vals)[ord][(length(ord) - n_group + 1):length(ord)]
  random_ids <- sample(setdiff(names(vals), c(low_ids, high_ids)), n_group)

  # Output file
  out_file <- file.path(outdir, paste0(trait, ".txt"))
  write.table(
    data.frame(
      snpid=character(), chr=character(), pos=integer(),
      high_ref=integer(), high_alt=integer(),
      low_ref=integer(), low_alt=integer(),
      random_ref=integer(), random_alt=integer()
    ),
    out_file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE
  )

  ### SNP chunk loop
  for (i in seq(1, n_snps, by = chunk_size)) {

    idx <- i:min(i + chunk_size - 1, n_snps)
    snp_chunk <- snp_ids[idx]

    dt <- fread(
      geno_file,
      select = c("ID", snp_chunk)
    )

    mat <- as.matrix(dt[, -1, with = FALSE])
    rownames(mat) <- dt$ID

    sub <- function(ids) mat[ids, , drop=FALSE]

    count <- function(m) {
      ref <- colSums(m == 0) * 2 + colSums(m == 1)
      alt <- colSums(m == 2) * 2 + colSums(m == 1)
      list(ref=ref, alt=alt)
    }

    h <- count(sub(high_ids))
    l <- count(sub(low_ids))
    r <- count(sub(random_ids))

    # Parse chr/pos
    chr <- sub("_.*", "", snp_chunk)
    pos <- as.integer(sub(".*_", "", snp_chunk))

    out <- data.frame(
      snpid = snp_chunk,
      chr = chr,
      pos = pos,
      high_ref = h$ref,
      high_alt = h$alt,
      low_ref = l$ref,
      low_alt = l$alt,
      random_ref = r$ref,
      random_alt = r$alt
    )

    fwrite(out, out_file, sep="\t", append=TRUE)
  }
}
