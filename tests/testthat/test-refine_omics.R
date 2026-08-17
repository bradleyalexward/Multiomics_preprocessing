testthat::test_that("prepare functions expose impute_behaviour choices", {
  testthat::expect_true("impute_behaviour" %in% names(formals(refineMetabolomics)))
  testthat::expect_false("impute_by_batch" %in% names(formals(refineMetabolomics)))
  testthat::expect_identical(
    eval(formals(refineMetabolomics)$impute_behaviour),
    c("per batch", "global")
  )

  testthat::expect_true("impute_behaviour" %in% names(formals(refineProteomics)))
  testthat::expect_identical(
    eval(formals(refineProteomics)$impute_behaviour),
    c("global", "per batch")
  )
})

testthat::test_that("prepare functions validate impute_behaviour values early", {
  testthat::expect_error(
    refineMetabolomics(impute_behaviour = "bad choice"),
    "one of"
  )
  testthat::expect_error(
    refineProteomics(impute_behaviour = "bad choice"),
    "one of"
  )
})

testthat::test_that("refineProteomics can impute per batch with non-mixed methods", {
  testthat::skip_if_not_installed("QFeatures")

  sample_meta <- data.frame(
    Sample = c("S1", "S2", "S3", "S4"),
    Batch = c("B1", "B1", "B2", "B2"),
    stringsAsFactors = FALSE
  )

  precursors <- data.frame(
    Precursor.Id = c("pep1", "pep2"),
    Protein.Group = c("P1", "P2"),
    S1 = c(10, 5),
    S2 = c(NA, 6),
    S3 = c(11, 7),
    S4 = c(12, NA),
    stringsAsFactors = FALSE
  )

  proteins <- refineProteomics(
    input_files = list(precursors),
    sample_meta_data = sample_meta,
    precursor_abundance_columns = c("S1", "S2", "S3", "S4"),
    impute_method = "with",
    impute_with_value = 1,
    impute_behaviour = "per batch",
    normalisation_method = "center.mean",
    protein_aggregator_method = base::colMeans,
    return_mnar_results = FALSE,
    verbose = FALSE
  )

  # Aggregation was not requested, so the result is the always-list contract
  # with the abundance matrix under $abundances (and $precursors as an alias).
  testthat::expect_type(proteins, "list")
  testthat::expect_s3_class(proteins$abundances, "data.frame")
  testthat::expect_false(anyNA(proteins$abundances))
  testthat::expect_equal(sort(colnames(proteins$abundances)), c("S1", "S2", "S3", "S4"))
  testthat::expect_identical(proteins$abundances, proteins$precursors)
  testthat::expect_identical(proteins$final_assay_name, "normalised")
})

testthat::test_that("aggregation switches the result alias and assay name", {
  testthat::skip_if_not_installed("QFeatures")

  sample_meta <- data.frame(Sample = c("S1", "S2", "S3", "S4"),
                            stringsAsFactors = FALSE)
  precursors <- data.frame(
    Precursor.Id = c("pep1", "pep2", "pep3", "pep4"),
    Protein.Group = c("P1", "P1", "P2", "P2"),
    S1 = c(10, 5, 8, 9), S2 = c(11, 6, 7, 8),
    S3 = c(12, 7, 9, 10), S4 = c(9, 8, 6, 11),
    stringsAsFactors = FALSE
  )

  res <- refineProteomics(
    input_files = list(precursors),
    sample_meta_data = sample_meta,
    precursor_abundance_columns = c("S1", "S2", "S3", "S4"),
    impute_method = "none",
    normalisation_method = "center.mean",
    protein_aggregator_method = base::colMeans,
    protein_aggregator_column = "Protein.Group",
    return_mnar_results = FALSE,
    verbose = FALSE
  )

  testthat::expect_identical(res$final_assay_name, "aggregated")
  testthat::expect_identical(res$abundances, res$proteins)
  testthat::expect_null(res$precursors)
  testthat::expect_equal(sort(rownames(res$proteins)), c("P1", "P2"))
})

testthat::test_that("refineProteomics blocks per-batch mixed imputation when a batch is fully missing", {
  testthat::skip_if_not_installed("QFeatures")
  testthat::skip_if_not_installed("missForest")
  testthat::skip_if_not_installed("imputeLCMD")

  sample_meta <- data.frame(
    Sample = c("S1", "S2", "S3", "S4"),
    Batch = c("B1", "B1", "B2", "B2"),
    Condition = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  precursors <- data.frame(
    Precursor.Id = c("pep1", "pep2"),
    Protein.Group = c("P1", "P2"),
    S1 = c(10, NA),
    S2 = c(12, NA),
    S3 = c(11, 20),
    S4 = c(13, 22),
    stringsAsFactors = FALSE
  )

  testthat::expect_error(
    refineProteomics(
      input_files = list(precursors),
      sample_meta_data = sample_meta,
      precursor_abundance_columns = c("S1", "S2", "S3", "S4"),
      mnar_variables = "Condition",
      mnar_significance_threshold = 1,
      missing_feature_proportion_threshold = c(1, 1),
      impute_method = "mixed",
      impute_behaviour = "per batch",
      normalisation_method = "center.mean",
      protein_aggregator_method = base::colMeans,
      verbose = FALSE
    ),
    "impute_behaviour = 'global'"
  )
})

testthat::test_that("missing optional imputation engines are reported before any work", {
  testthat::expect_identical(
    OmicsRefinery:::mp_impute_engines("mixed"),
    c("missForest", "imputeLCMD")
  )
  testthat::expect_identical(
    OmicsRefinery:::mp_impute_engines("QRILC"), "imputeLCMD")
  testthat::expect_identical(
    OmicsRefinery:::mp_impute_engines("zero"), character(0))
})
