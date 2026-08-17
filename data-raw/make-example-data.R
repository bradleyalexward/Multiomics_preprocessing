# Generate the synthetic example data shipped in inst/extdata.
#
# THESE FILES CONTAIN NO REAL MEASUREMENTS AND NO HUMAN SUBJECT DATA.
# Every value is drawn from a random number generator by this script. The
# column names and file layout imitate the shape of a cleaned Progenesis QI
# peak table (metabolomics) and a wide precursor table exported from DIA-NN
# (proteomics), so that the worked examples in the vignette look like real
# input, but nothing here is derived from any measured sample.
#
# Run from the package root with:
#   Rscript data-raw/make-example-data.R

set.seed(20260817)

out_dir <- file.path("inst", "extdata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 30 samples rather than a token handful: mixed imputation fits a distribution
# per sample within each block, and a smaller design cannot demonstrate that
# honestly.
n_samples <- 30
samples   <- sprintf("S%02d", seq_len(n_samples))

# Balanced design: group and sex are crossed with batch, so batch correction
# has something to remove without being confounded with the biology.
sample_meta <- data.frame(
  Sample_ID         = samples,
  Sample_group      = rep(rep(c("A", "B"), each = 5), 3),
  Sex               = rep(c("M", "F"), n_samples / 2),
  Metabolomic_batch = rep(c("Batch 1", "Batch 2"), each = 15),
  Proteomic_batch   = rep(c("Batch 1", "Batch 2"), each = 15),
  stringsAsFactors  = FALSE
)

group_is_b <- sample_meta$Sample_group == "B"

# ---------------------------------------------------------------------------
# Shared feature simulator
# ---------------------------------------------------------------------------
# Returns a features x samples matrix of positive abundances on the linear
# scale, with three kinds of missingness:
#   "complete"  - detected nearly everywhere
#   "mar"       - missing at random, unrelated to any covariate
#   "mnar"      - detection depends on Sample_group (what mixed imputation
#                 is meant to recognise)
simulate_block <- function(n_features, batch, prefix) {
  kind <- sample(c("complete", "mar", "mnar"),
                 n_features, replace = TRUE, prob = c(0.60, 0.15, 0.25))

  baseline <- stats::rnorm(n_features, mean = 20, sd = 2)          # log2 scale
  batch_shift <- stats::rnorm(n_features, mean = 0, sd = 0.9)      # batch effect
  group_shift <- ifelse(stats::runif(n_features) < 0.25,
                        stats::rnorm(n_features, 0, 1.2), 0)       # real biology

  x <- matrix(NA_real_, nrow = n_features, ncol = n_samples,
              dimnames = list(sprintf("%s_%04d", prefix, seq_len(n_features)), samples))

  for (i in seq_len(n_features)) {
    mu <- baseline[i] +
      batch_shift[i] * (batch == levels(factor(batch))[2]) +
      group_shift[i] * group_is_b
    x[i, ] <- 2^(mu + stats::rnorm(n_samples, 0, 0.35))

    # The MNAR contrast is strong enough to be detected reliably at the default
    # significance threshold, but not so extreme that group B has nothing left
    # for QRILC to estimate a per-sample distribution from. Real studies have
    # thousands of features and do not face that second constraint.
    detect_p <- switch(kind[i],
                       complete = rep(0.97, n_samples),
                       mar      = rep(0.75, n_samples),
                       mnar     = ifelse(group_is_b, 0.35, 0.95))
    keep <- stats::runif(n_samples) < detect_p

    # Batch-wise imputation needs at least one observed value per batch, so
    # force a minimum of two observations in each batch.
    for (b in unique(batch)) {
      in_b <- batch == b
      if (sum(keep[in_b]) < 2) {
        keep[which(in_b)[order(-x[i, in_b])[1:2]]] <- TRUE
      }
    }
    x[i, !keep] <- NA_real_
  }
  x
}

# ---------------------------------------------------------------------------
# Metabolomics: four peak tables (two ionisation modes x two runs)
# ---------------------------------------------------------------------------
elements <- function(n) {
  sprintf("C%dH%dN%dO%d",
          sample(6:30, n, TRUE), sample(6:50, n, TRUE),
          sample(0:5, n, TRUE),  sample(1:9, n, TRUE))
}

write_metabolite_file <- function(dataset, n_peaks, file) {
  x <- simulate_block(n_peaks, sample_meta$Metabolomic_batch, dataset)

  # Several peaks map to the same compound, so aggregation has work to do.
  n_compounds <- ceiling(n_peaks / 2.5)
  compound_id <- sprintf("%s_CMP_%03d", dataset,
                         sample(seq_len(n_compounds), n_peaks, replace = TRUE))

  df <- data.frame(
    Compound           = rownames(x),
    Compound_ID        = compound_id,
    Formula            = elements(n_peaks),
    Description        = paste("Synthetic feature", seq_len(n_peaks)),
    Retention_time_min = round(stats::runif(n_peaks, 0.5, 18), 3),
    mz                 = round(stats::runif(n_peaks, 90, 1000), 4),
    stringsAsFactors   = FALSE
  )

  # Peak tables report a zero for "not detected", not an empty cell.
  quant <- round(x, 1)
  quant[is.na(quant)] <- 0
  df <- cbind(df, as.data.frame(quant))

  utils::write.csv(df, file, row.names = FALSE)
  invisible(df)
}

write_metabolite_file("NegM1", 200, file.path(out_dir, "Metabolites_NegM1.csv"))
write_metabolite_file("NegM2", 200, file.path(out_dir, "Metabolites_NegM2.csv"))
write_metabolite_file("PosM1", 200, file.path(out_dir, "Metabolites_PosM1.csv"))
write_metabolite_file("PosM2", 200, file.path(out_dir, "Metabolites_PosM2.csv"))

# ---------------------------------------------------------------------------
# Proteomics: three fractions in the wide layout prepare_proteins() expects
# ---------------------------------------------------------------------------
# Note the deliberate schema differences: each fraction carries one annotation
# column the others do not. prepare_proteins() combines these by taking the
# union of annotation columns.
write_protein_file <- function(fraction, n_precursors, extra, file) {
  x <- simulate_block(n_precursors, sample_meta$Proteomic_batch, fraction)

  n_groups <- ceiling(n_precursors / 3)
  group_no <- sample(seq_len(n_groups), n_precursors, replace = TRUE)

  df <- data.frame(
    Precursor.Id     = sprintf("%s_PEP%04d_2", fraction, seq_len(n_precursors)),
    Protein.Group    = sprintf("SYN%04d", group_no),
    Genes            = sprintf("GENE%04d", group_no),
    stringsAsFactors = FALSE
  )
  df[[names(extra)]] <- extra[[1]](n_precursors)

  # Search engines leave undetected precursors empty rather than zero.
  df <- cbind(df, as.data.frame(round(x, 1)))

  utils::write.csv(df, file, row.names = FALSE)
  invisible(df)
}

write_protein_file(
  "Frac1", 400,
  list(Modified.Sequence = function(n) sprintf("(UniMod:4)PEPTIDE%03dK", seq_len(n))),
  file.path(out_dir, "Proteomics_fraction_1.csv")
)
write_protein_file(
  "Frac2", 400,
  list(Proteotypic = function(n) sample(c(0L, 1L), n, TRUE)),
  file.path(out_dir, "Proteomics_fraction_2.csv")
)
write_protein_file(
  "Frac3", 400,
  list(Precursor.Charge = function(n) sample(2:4, n, TRUE)),
  file.path(out_dir, "Proteomics_fraction_3.csv")
)

# ---------------------------------------------------------------------------
utils::write.csv(sample_meta, file.path(out_dir, "samples.csv"), row.names = FALSE)

message("Wrote ", length(list.files(out_dir)), " synthetic example files to ", out_dir)
