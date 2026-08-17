# Regression tests for defects that produced no error message.
#
# Each test below pins a bug that changed results silently. Failures that
# announce themselves (bad column names, duplicate feature IDs, logs of zero)
# are deliberately not duplicated here -- if those regress, the code stops.

samples_12 <- paste0("S", sprintf("%02d", 1:12))

make_meta <- function(sample_col = "Sample") {
  m <- data.frame(
    ID        = samples_12,
    Batch     = rep(c("B1", "B2"), each = 6),
    Condition = rep(c("A", "B"), each = 6),
    Sex       = rep(c("M", "F"), 6),
    stringsAsFactors = FALSE
  )
  names(m)[1] <- sample_col
  m
}

make_precursors <- function(n = 8, seed = 1) {
  set.seed(seed)
  m <- matrix(round(stats::runif(n * length(samples_12), 500, 5000), 1),
              nrow = n, dimnames = list(NULL, samples_12))
  data.frame(
    Precursor.Id  = sprintf("pep%02d", seq_len(n)),
    Protein.Group = sprintf("P%02d", rep(seq_len(n / 2), each = 2)),
    m, check.names = FALSE, stringsAsFactors = FALSE
  )
}

# --- ordered column selection ----------------------------------------------

testthat::test_that("abundance columns are built in the requested order, not file order", {
  testthat::skip_if_not_installed("QFeatures")

  reversed <- rev(samples_12)

  res <- refineMetabolomics(
    input_files = list(stats::setNames(make_precursors(), c("Compound", "Compound_ID", samples_12))),
    peak_abundance_columns = reversed,
    peak_id_column = "Compound",
    sample_meta_data = make_meta(),
    impute_method = "none",
    normalisation_method = "center.mean",
    return_mnar_results = FALSE,
    verbose = FALSE
  )

  # A logical mask would silently hand back file order here.
  testthat::expect_identical(
    colnames(SummarizedExperiment::assay(res$qfeatures[["raw"]])),
    reversed
  )
  testthat::expect_identical(colnames(res$abundances), reversed)
})

testthat::test_that("shuffled metadata rows give identical results to ordered metadata", {
  testthat::skip_if_not_installed("QFeatures")

  args <- list(
    input_files = list(make_precursors()),
    precursor_abundance_columns = samples_12,
    impute_method = "none",
    normalisation_method = "center.mean",
    return_mnar_results = FALSE,
    verbose = FALSE
  )

  ordered  <- do.call(refineProteomics, c(args, list(sample_meta_data = make_meta())))
  set.seed(11)
  shuffled <- do.call(refineProteomics,
                      c(args, list(sample_meta_data = make_meta()[sample(12), ])))

  testthat::expect_equal(ordered$abundances, shuffled$abundances)
})

# --- missingness classification --------------------------------------------

classify <- function(detected, covariates = c("Condition", "Sex"), cutoff = 0.05) {
  df <- data.frame(Feature = rownames(detected),
                   ifelse(detected, 1000, NA_real_),
                   check.names = FALSE, stringsAsFactors = FALSE)
  se <- QFeatures::readSummarizedExperiment(df, quantCols = samples_12,
                                            fnames = "Feature")
  qf <- QFeatures::QFeatures(list(raw = se))
  info <- OmicsRefinery:::mp_build_sample_info(
    make_meta(), samples_12, "meta", "Sample")
  OmicsRefinery:::mp_classify_missingness(
    qf, info, covariates, "raw", "Feature", cutoff, 0.3)$results
}

testthat::test_that("missingness unrelated to the covariates is rarely called MNAR", {
  testthat::skip_if_not_installed("QFeatures")

  set.seed(3)
  n <- 200
  detected <- matrix(stats::runif(n * 12) < 0.8, nrow = n,
                     dimnames = list(sprintf("f%03d", seq_len(n)), samples_12))
  detected[rowSums(detected) == 0, ] <- TRUE   # keep every feature testable

  # Benjamini-Hochberg controls the false discovery rate, so the assertion is
  # that the rate stays near the nominal level -- not that it is exactly zero,
  # which no honest hypothesis test can promise.
  res <- classify(detected, cutoff = 0.01)
  testthat::expect_lt(mean(res$MNAR), 0.02)
})

testthat::test_that("a perfectly separated feature is called MNAR", {
  testthat::skip_if_not_installed("QFeatures")

  # Detected in every Condition A sample and no Condition B sample: the
  # strongest possible signal. The old Wald-based test scored this p ~ 1.
  detected <- matrix(TRUE, nrow = 2, ncol = 12,
                     dimnames = list(c("sep", "flat"), samples_12))
  detected["sep", 7:12] <- FALSE

  res <- classify(detected, covariates = "Condition")
  testthat::expect_true(res$MNAR[res$Feature == "sep"])

  # A feature detected everywhere carries no signal and must stay MAR.
  testthat::expect_false(res$MNAR[res$Feature == "flat"])
  testthat::expect_equal(res$pvals[res$Feature == "flat"], 1)
})

testthat::test_that("one multiple-testing correction is applied, across features", {
  testthat::skip_if_not_installed("QFeatures")

  detected <- matrix(TRUE, nrow = 3, ncol = 12,
                     dimnames = list(c("a", "b", "c"), samples_12))
  detected["a", 7:12] <- FALSE
  detected["b", c(1, 8)] <- FALSE
  detected["c", 10:12] <- FALSE

  res <- classify(detected, covariates = "Condition", cutoff = 0.05)
  testthat::expect_equal(res$adj.pvals,
                         stats::p.adjust(res$pvals, method = "BH"))
})

# --- abundance coercion ----------------------------------------------------

testthat::test_that("factor abundance columns keep their displayed values", {
  df <- data.frame(id = c("a", "b", "c"),
                   S1 = factor(c("10", "100", "20")),
                   stringsAsFactors = FALSE)
  out <- OmicsRefinery:::mp_coerce_abundance_to_numeric(df, 2L, "test")

  # as.numeric() on a factor returns level codes: this used to yield 1, 2, 3.
  testthat::expect_identical(out$S1, c(10, 100, 20))

  bad <- data.frame(id = "a", S1 = factor("not a number"), stringsAsFactors = FALSE)
  testthat::expect_error(
    OmicsRefinery:::mp_coerce_abundance_to_numeric(bad, 2L, "test"),
    "Non-numeric values"
  )
})

# --- sample identity -------------------------------------------------------

testthat::test_that("ambiguous sample ID columns error, and sample_id_column resolves it", {
  meta <- make_meta()
  meta$Paired_ID <- rev(meta$Sample)          # a second full-coverage column

  testthat::expect_error(
    OmicsRefinery:::mp_build_sample_info(meta, samples_12, "meta"),
    "Multiple columns"
  )

  info <- OmicsRefinery:::mp_build_sample_info(
    meta, samples_12, "meta", sample_id_column = "Sample")
  testthat::expect_identical(info$Sample, samples_12)
  testthat::expect_identical(attr(info, "sample_col"), "Sample")
})

testthat::test_that("tibble metadata survives and keeps the QC plots", {
  testthat::skip_if_not_installed("QFeatures")
  testthat::skip_if_not_installed("tibble")

  # recordPlot() needs a device with the display list enabled.
  grDevices::pdf(NULL)
  grDevices::dev.control(displaylist = "enable")
  on.exit(grDevices::dev.off(), add = TRUE)

  res <- refineProteomics(
    input_files = list(make_precursors(n = 20)),
    sample_meta_data = tibble::as_tibble(make_meta()),
    precursor_abundance_columns = samples_12,
    impute_method = "none",
    normalisation_method = "center.mean",
    batch_correction = TRUE,
    batch_column = "Batch",
    plot_qc_charts = TRUE,
    return_mnar_results = FALSE,
    verbose = FALSE
  )

  # Row names do not survive dplyr verbs on a tibble, which used to make both
  # PCA plots disappear with only a log line to say so.
  testthat::expect_true(!is.null(res$plots$before_batch_correction))
  testthat::expect_true(!is.null(res$plots$after_batch_correction))
  testthat::expect_s3_class(res$plots$before_batch_correction, "ggplot")
})

# --- thresholds ------------------------------------------------------------

testthat::test_that("a scalar threshold is never stricter for MNAR than for MAR", {
  mnar_thr <- OmicsRefinery:::mp_mnar_threshold

  testthat::expect_equal(mnar_thr(0.3), 0.75)   # documented allowance
  testthat::expect_equal(mnar_thr(0.9), 0.9)    # never below the MAR threshold
  testthat::expect_equal(mnar_thr(c(0.3, 0.6)), 0.6)
})

testthat::test_that("features exactly on the threshold are retained in both paths", {
  testthat::skip_if_not_installed("QFeatures")

  # 3 of 12 samples missing is exactly 0.25; the threshold is set to 0.25 so
  # the feature sits on the boundary. "<=" keeps it in both branches.
  df <- make_precursors(n = 4)
  df[1, samples_12[1:3]] <- NA

  res <- refineProteomics(
    input_files = list(df),
    sample_meta_data = make_meta(),
    precursor_abundance_columns = samples_12,
    missing_feature_proportion_threshold = 0.25,
    impute_method = "none",
    normalisation_method = "center.mean",
    return_mnar_results = FALSE,
    verbose = FALSE
  )
  testthat::expect_true("pep01" %in% rownames(res$abundances))
})

# --- session hygiene -------------------------------------------------------

testthat::test_that("the caller's random number generator is restored", {
  testthat::skip_if_not_installed("QFeatures")

  # Build the fixture first: make_precursors() seeds the RNG itself, and as a
  # lazily evaluated argument it would run inside the call being tested.
  fixture <- make_precursors()
  meta <- make_meta()

  set.seed(99)
  before_kind <- RNGkind()
  before_seed <- get(".Random.seed", envir = globalenv())

  invisible(refineProteomics(
    input_files = list(fixture),
    sample_meta_data = meta,
    precursor_abundance_columns = samples_12,
    impute_method = "with",
    impute_with_value = 1,
    imputation_seeds = 1L,
    normalisation_method = "center.mean",
    return_mnar_results = FALSE,
    verbose = FALSE
  ))

  testthat::expect_identical(RNGkind(), before_kind)
  testthat::expect_identical(get(".Random.seed", envir = globalenv()), before_seed)
})

# --- fraction combination --------------------------------------------------

testthat::test_that("fractions with different annotation columns combine", {
  testthat::skip_if_not_installed("QFeatures")

  f1 <- make_precursors(n = 6, seed = 1)
  f2 <- make_precursors(n = 6, seed = 2)
  f1$Modified.Sequence <- paste0("MOD", seq_len(6))   # only in fraction 1
  f2$Proteotypic <- rep(1L, 6)                        # only in fraction 2

  res <- refineProteomics(
    input_files = list(f1, f2),
    fraction_file_names = c("F1", "F2"),
    sample_meta_data = make_meta(),
    precursor_abundance_columns = samples_12,
    impute_method = "none",
    normalisation_method = "center.mean",
    return_mnar_results = FALSE,
    verbose = FALSE
  )

  # base rbind() would have stopped with "names do not match previous names".
  testthat::expect_equal(nrow(res$abundances), 12)
  rd <- SummarizedExperiment::rowData(res$qfeatures[["raw"]])
  testthat::expect_true(all(c("Modified.Sequence", "Proteotypic") %in% colnames(rd)))
  testthat::expect_equal(sum(is.na(rd$Proteotypic)), 6)
})

# --- numerical safety ------------------------------------------------------

testthat::test_that("a log transform of non-positive values stops", {
  testthat::skip_if_not_installed("QFeatures")

  df <- make_precursors(n = 4)
  df[1, samples_12[1]] <- 0

  testthat::expect_error(
    refineProteomics(
      input_files = list(df),
      sample_meta_data = make_meta(),
      precursor_abundance_columns = samples_12,
      zero_handling = "leave",
      impute_method = "none",
      normalisation_method = c("log2", "center.mean"),
      return_mnar_results = FALSE,
      verbose = FALSE
    ),
    "zero or negative"
  )
})

# --- end-to-end ------------------------------------------------------------

testthat::test_that("the bundled metabolomics example runs end to end", {
  testthat::skip_if_not_installed("QFeatures")
  testthat::skip_if_not_installed("missForest")
  testthat::skip_if_not_installed("imputeLCMD")
  testthat::skip_on_cran()

  dir <- system.file("extdata", package = "OmicsRefinery")
  testthat::skip_if(dir == "", "example data not installed")

  samples <- utils::read.csv(file.path(dir, "samples.csv"), check.names = FALSE)
  neg <- utils::read.csv(file.path(dir, "Metabolites_NegM1.csv"), check.names = FALSE)
  pos <- utils::read.csv(file.path(dir, "Metabolites_PosM1.csv"), check.names = FALSE)
  cols <- intersect(samples$Sample_ID, colnames(neg))

  res <- refineMetabolomics(
    input_files = list(neg, pos),
    dataset_file_names = c("NegM1", "PosM1"),
    peak_abundance_columns = cols,
    peak_id_column = "Compound",
    sample_meta_data = samples,
    impute_method = "mixed",
    mnar_variables = c("Sample_group", "Sex"),
    imputation_seeds = 1L,
    impute_behaviour = "per batch",
    normalisation_method = "vsn",
    batch_correction = TRUE,
    batch_column = "Metabolomic_batch",
    batch_correction_variable_columns = c("Sample_group", "Sex"),
    metabolite_aggregator_column = "Compound_ID",
    verbose = FALSE
  )

  testthat::expect_identical(res$final_assay_name, "aggregated")
  testthat::expect_identical(colnames(res$metabolites), cols)
  testthat::expect_true(all(is.finite(as.matrix(res$metabolites))))
  testthat::expect_true(nrow(res$mnar_results[["raw"]]) > 0)
})

testthat::test_that("the bundled proteomics example runs end to end", {
  testthat::skip_if_not_installed("QFeatures")
  testthat::skip_if_not_installed("missForest")
  testthat::skip_if_not_installed("imputeLCMD")
  testthat::skip_on_cran()

  dir <- system.file("extdata", package = "OmicsRefinery")
  testthat::skip_if(dir == "", "example data not installed")

  samples <- utils::read.csv(file.path(dir, "samples.csv"), check.names = FALSE)
  fr <- lapply(1:3, function(i) {
    utils::read.csv(file.path(dir, sprintf("Proteomics_fraction_%d.csv", i)),
                    check.names = FALSE)
  })
  cols <- intersect(samples$Sample_ID, colnames(fr[[1]]))

  res <- refineProteomics(
    input_files = fr,
    fraction_file_names = c("F1", "F2", "F3"),
    precursor_abundance_columns = cols,
    sample_meta_data = samples,
    impute_method = "mixed",
    mnar_variables = c("Sample_group", "Sex"),
    imputation_seeds = 2L,
    impute_behaviour = "global",
    batch_specific_features_threshold = 0.9,
    normalisation_method = "vsn",
    batch_correction = TRUE,
    batch_column = "Proteomic_batch",
    batch_correction_variable_columns = c("Sample_group", "Sex"),
    protein_aggregator_column = "Protein.Group",
    verbose = FALSE
  )

  testthat::expect_identical(res$final_assay_name, "aggregated")
  testthat::expect_identical(colnames(res$proteins), cols)
  testthat::expect_true(all(is.finite(as.matrix(res$proteins))))
  testthat::expect_length(res$mnar_results, 3)
})
