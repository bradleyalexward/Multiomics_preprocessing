testthat::test_that("prepare functions expose impute_behaviour choices", {
  testthat::expect_true("impute_behaviour" %in% names(formals(prepare_metabolites)))
  testthat::expect_false("impute_by_batch" %in% names(formals(prepare_metabolites)))
  testthat::expect_identical(
    eval(formals(prepare_metabolites)$impute_behaviour),
    c("per batch", "global")
  )

  testthat::expect_true("impute_behaviour" %in% names(formals(prepare_proteins)))
  testthat::expect_identical(
    eval(formals(prepare_proteins)$impute_behaviour),
    c("global", "per batch")
  )
})

testthat::test_that("prepare functions validate impute_behaviour values early", {
  testthat::expect_error(
    prepare_metabolites(impute_behaviour = "bad choice"),
    "one of"
  )
  testthat::expect_error(
    prepare_proteins(impute_behaviour = "bad choice"),
    "one of"
  )
})

testthat::test_that("prepare_proteins can impute per batch with non-mixed methods", {
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("MsCoreUtils")
  testthat::skip_if_not_installed("QFeatures")
  testthat::skip_if_not_installed("SummarizedExperiment")

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

  proteins <- prepare_proteins(
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

  testthat::expect_s3_class(proteins, "data.frame")
  testthat::expect_false(anyNA(proteins))
  testthat::expect_equal(sort(colnames(proteins)), c("S1", "S2", "S3", "S4"))
})

testthat::test_that("prepare_proteins blocks per-batch mixed imputation when a batch is fully missing", {
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("MsCoreUtils")
  testthat::skip_if_not_installed("QFeatures")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("tibble")
  testthat::skip_if_not_installed("tidyr")

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
    prepare_proteins(
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
