#' Prepare proteomics data from precursor intensities to a protein abundance matrix
#'
#' @param input_files List of one or more precursor-level data.frames.
#' Example: list(fraction_1, fraction_2)
#' @param fraction_file_names Optional assay names, one per input file.
#' Example: c("fraction_1", "fraction_2")
#' @param precursor_abundance_columns Abundance columns (numeric indices, names, or logical selector).
#' Example: 5:70 or c("Quant1", "Quant2", "Quant3", etc)
#' @param precursor_id_column Feature identifier column in each input table (must share same name).
#' @param sample_meta_data Sample metadata data.frame (one row per sample).
#' @param sample_id_column Optional column in `sample_meta_data` holding the sample
#' identifiers that match the abundance column names. If `NULL` (the default) the
#' column is detected automatically; detection requires exactly one column to
#' contain every assay sample ID and errors if several do.
#' @param zero_handling Either `"na"` (replace 0 with NA) or `"leave"`.
#' Default is "leave"
#' @param impute_method Missing-value imputation strategy.
#' Can be any of: "none", "mixed", "bpca", "knn", "QRILC", "MLE", "MLE2", "MinDet", "MinProb", "min", "zero", "nbavg", "with", "RF".
#' When mixed is selected, MNAR features imputed via QRILC, MAR features imputed via RF.
#' Several methods rely on optional packages (see Details); these are checked before
#' any work is done and reported with an install command if missing.
#' @param mnar_variables Optional character vector of column name(s) in `sample_meta_data`
#' used when `impute_method = "mixed"`.
#' These covariate columns are used for logistic regression-based MNAR/MAR determination and
#' must not include the sample ID column.
#' Example: `"Sex"` or `c("Sex", "Age")`
#' @param mnar_significance_threshold MNAR adjusted p-value threshold for mixed imputation.
#' Default 0.01
#' @param return_mnar_results Return MNAR classification tables.
#' TRUE FALSE whether to return data table with statistical results of MNAR logistic regression analysis.
#' Default is TRUE
#' @param impute_with_value Replacement value used only when `impute_method = "with"`.
#' Must be numeric.
#' @param imputation_seeds Optional integer seed or integer vector controlling stochastic
#' imputation steps. `NULL` keeps the current RNG behaviour.
#' If length 1, successive imputation steps use `seed`, `seed + 1`, `seed + 2`, etc.
#' If length > 1, it must match the number of imputation steps in the current run.
#' The caller's random number generator state is restored when the function exits.
#' @param impute_behaviour Whether to impute missing values globally across all samples (`"global"`)
#' or separately within each batch using `batch_column` (`"per batch"`).
#' Default is `"global"`. For `impute_method = "mixed"`, retained features must have at least one observed value in each batch when `"per batch"` is used.
#' @param missing_feature_proportion_threshold Length 1 or 2 numeric threshold(s) in `[0, 1]`.
#' Default is c(0.3, 0.75). Here, features with more than 30% missingness are removed (first numeric) unless they are MNAR, in which case they are retained up to 75% missingness.
#' If only one numeric is supplied like c(0.3), it sets the threshold for MAR features;
#' MNAR features are then retained up to 0.75, or up to the supplied value when that is
#' higher than 0.75. Features are retained when missingness is less than or equal to the
#' applicable threshold.
#' Two numbers must only be supplied if impute_method = "mixed"
#' @param batch_specific_features_threshold Optional missingness threshold used to remove features mostly missing within any batch before imputation.
#' Default is `NULL`, meaning no batch-specific feature removal is performed.
#' If a number is supplied, features with more than this proportion of missing values in any batch are removed.
#' Applied to each fraction separately, before fractions are combined.
#' Example: `0.9` removes features with more than 90% missingness in at least one batch.
#' @param normalisation_method Length 1 normalisation method, or length 2 transform+method. Transform methods include log, log2, log10.
#' Normalise methods include "sum", "max", "center.mean", "center.median", "div.mean", "div.median", "diff.median", "quantiles", "quantiles.robust" or "vsn".
#' Example: c("vsn") or c("log10", "center.median").
#' Default is vsn.
#' @param batch_correction Apply ComBat batch correction.
#' @param batch_column Column in `sample_meta_data` identifying which sample belongs to which batch.
#' Used for optional batch-specific filtering, optional batch-wise imputation, and ComBat batch correction.
#' Must be supplied explicitly when `batch_specific_features_threshold` is used.
#' @param batch_correction_variable_columns Optional biological covariates of interest to preserve during batch correction
#' @param protein_aggregator_method Protein aggregation function. Function by which peptides are aggregated to protein level when aggregation is requested.
#' Default is MsCoreUtils::robustSummary, other examples include MsCoreUtils::medianPolish, base::colMeans, matrixStats::colMedians, etc.
#' Function must take a matrix as input and return a vector of length equal to the column length of the input.
#' It must also accept an `na.rm` argument (or `...`), because aggregation is always called with `na.rm = TRUE`.
#' @param protein_aggregator_column Optional column in precursor-level data.frames mapping each precursor to their identified protein (must share same name).
#' If `NULL`, aggregation is skipped.
#' @param plot_qc_charts Return QC plots including density plots of raw, normalised, batch-corrected, and aggregated abundances when available, plus PCA before and after batch correction.
#' QC plots are captured with [grDevices::recordPlot()] and are intended for interactive
#' sessions; on a file device they may return blank.
#' @param show_plot_legends Legend behavior for density plots (`"auto"`, `"show"`, `"hide"`).
#' Default is auto, which hides legends if more than 12 samples are present.
#' @param verbose If `TRUE`, print step-by-step progress messages.
#' Default is TRUE.
#'
#' @details
#' # Input layout
#'
#' Each element of `input_files` must be a *wide* table: one row per precursor and
#' one column per sample. Long-format search engine output (for example a DIA-NN
#' `report.tsv`) must be pivoted to that shape first.
#'
#' Fractions may carry different ancillary annotation columns; they are combined by
#' taking the union, with `NA` where a fraction lacks a column. Only
#' `precursor_id_column`, `protein_aggregator_column` and the abundance columns must
#' be common to every fraction.
#'
#' # Missingness classification
#'
#' When `impute_method = "mixed"`, each precursor is tested with a logistic regression
#' asking whether its *detection* pattern is explained by `mnar_variables`. The test
#' is a likelihood-ratio comparison of the covariate model against an intercept-only
#' model, with a single Benjamini-Hochberg correction across features.
#'
#' This detects association between detection and the metadata you supplied. It is a
#' practical heuristic for routing features to left-censored (QRILC) imputation, not a
#' proof of the formal missing-not-at-random mechanism, which depends on values that
#' were never observed.
#'
#' # Optional dependencies
#'
#' Some imputation methods are implemented by packages that `MsCoreUtils` only
#' suggests, so they are not installed automatically:
#' `"mixed"` needs missForest and imputeLCMD; `"RF"` needs missForest;
#' `"QRILC"`, `"MinDet"` and `"MinProb"` need imputeLCMD; `"knn"` needs impute;
#' `"bpca"` needs pcaMethods; `"MLE"` needs norm. Quantile normalisation needs
#' preprocessCore.
#'
#' @return A list containing the final abundance matrix, the final assay name, the full
#' QFeatures workflow object, and optional extras such as MNAR tables and QC plots.
#' @export
refineProteomics <- function(input_files,
                             fraction_file_names = NULL,
                             precursor_abundance_columns,
                             precursor_id_column = "Precursor.Id",
                             sample_meta_data,
                             sample_id_column = NULL,
                             zero_handling = c("leave", "na"),
                             impute_method = c("none", "mixed", "bpca", "knn", "QRILC", "MLE", "MLE2", "MinDet", "MinProb", "min", "zero", "nbavg", "with", "RF"),
                             mnar_variables = NULL,
                             mnar_significance_threshold = 0.01,
                             return_mnar_results = TRUE,
                             impute_with_value = NULL,
                             imputation_seeds = NULL,
                             impute_behaviour = c("global", "per batch"),
                             missing_feature_proportion_threshold = c(0.3, 0.75),
                             batch_specific_features_threshold = NULL,
                             normalisation_method = c("vsn"),
                             batch_correction = FALSE,
                             batch_column = "Batch",
                             batch_correction_variable_columns = NULL,
                             protein_aggregator_method = MsCoreUtils::robustSummary,
                             protein_aggregator_column = NULL,
                             plot_qc_charts = FALSE,
                             show_plot_legends = c("auto", "show", "hide"),
                             verbose = TRUE) {

  batch_column_missing <- missing(batch_column)

  zero_handling <- match.arg(zero_handling)
  show_plot_legends <- match.arg(show_plot_legends)
  impute_method <- match.arg(impute_method)
  impute_behaviour <- match.arg(impute_behaviour)

  log_step <- function(fmt, ...) mp_log(verbose, fmt, ...)

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  mp_validate_common_args(input_files, sample_meta_data, verbose, batch_correction,
                          plot_qc_charts, return_mnar_results,
                          mnar_significance_threshold, impute_method,
                          impute_with_value)
  mp_assert_string(precursor_id_column, "precursor_id_column")
  mp_assert_string(sample_id_column, "sample_id_column", allow_null = TRUE)
  mp_validate_aggregator(protein_aggregator_column, protein_aggregator_method,
                         "protein_aggregator_column", "protein_aggregator_method")
  mp_validate_dataset_names(fraction_file_names, length(input_files),
                            "fraction_file_names")
  mp_validate_thresholds(missing_feature_proportion_threshold,
                         batch_specific_features_threshold, batch_column_missing)
  mp_validate_normalisation(normalisation_method)
  mp_validate_mnar_variables(mnar_variables, impute_method)
  imputation_seeds <- mp_validate_seeds(imputation_seeds)

  # Imputation switches the RNG kind; put the caller's session back afterwards.
  old_rng <- mp_save_rng()
  on.exit(mp_restore_rng(old_rng), add = TRUE)

  use_batchwise_imputation <- (impute_method != "none") && identical(impute_behaviour, "per batch")
  needs_batch_info <- !is.null(batch_specific_features_threshold) || use_batchwise_imputation || isTRUE(batch_correction)
  mp_validate_batch_args(needs_batch_info, batch_column, batch_correction,
                         batch_correction_variable_columns, sample_meta_data)

  # -------------------------------------------------------------------------
  # Build and clean inputs
  # -------------------------------------------------------------------------
  dfs <- input_files
  if (is.null(fraction_file_names)) {
    fraction_file_names <- if (length(dfs) == 1) "se" else paste0("Fraction", seq_along(dfs))
  }

  quant_cols_list <- vector("list", length(dfs))
  for (i in seq_along(dfs)) {
    frac <- fraction_file_names[i]
    if (!(precursor_id_column %in% colnames(dfs[[i]]))) {
      stop("`precursor_id_column` ('", precursor_id_column, "') not found in fraction '", frac, "'.")
    }
    quant_cols_list[[i]] <- mp_resolve_quant_cols(dfs[[i]], precursor_abundance_columns,
                                                  frac, "precursor_abundance_columns")
    protected_cols <- colnames(dfs[[i]])[quant_cols_list[[i]]]
    if (precursor_id_column %in% protected_cols) {
      stop("`precursor_id_column` cannot be included in `precursor_abundance_columns` for fraction '", frac, "'.")
    }
    if (!is.null(protein_aggregator_column)) {
      if (!(protein_aggregator_column %in% colnames(dfs[[i]]))) {
        stop("`protein_aggregator_column` ('", protein_aggregator_column, "') not found in fraction '", frac, "'.")
      }
      if (protein_aggregator_column %in% protected_cols) {
        stop("`protein_aggregator_column` cannot be included in `precursor_abundance_columns` for fraction '", frac, "'.")
      }
    }
    # Check the identifiers the user actually supplied, before any fraction
    # suffix is added: the suffix only disambiguates across fractions, never
    # within one, and "NA_Fraction_1" is neither missing nor empty.
    mp_validate_feature_ids(dfs[[i]][[precursor_id_column]], precursor_id_column, frac)
    dfs[[i]] <- mp_coerce_abundance_to_numeric(dfs[[i]], quant_cols_list[[i]], frac)
  }

  sample_names <- colnames(dfs[[1]])[quant_cols_list[[1]]]
  if (anyDuplicated(sample_names) > 0) {
    stop("Abundance column names in the first fraction must be unique so other fractions can be aligned to the same order.")
  }

  if (length(dfs) > 1) {
    for (i in 2:length(dfs)) {
      this_samples <- colnames(dfs[[i]])[quant_cols_list[[i]]]
      if (anyDuplicated(this_samples) > 0) {
        stop("Abundance column names in fraction '", fraction_file_names[i],
             "' must be unique so columns can be aligned to the first fraction.")
      }

      missing_samples <- setdiff(sample_names, this_samples)
      extra_samples <- setdiff(this_samples, sample_names)
      if (length(missing_samples) > 0 || length(extra_samples) > 0) {
        stop("Abundance columns in fraction '", fraction_file_names[i],
             "' do not match the first fraction. Missing: ",
             if (length(missing_samples) > 0) paste(missing_samples, collapse = ", ") else "none",
             ". Extra: ",
             if (length(extra_samples) > 0) paste(extra_samples, collapse = ", ") else "none",
             ".")
      }

      if (!identical(this_samples, sample_names)) {
        aligned_quant_cols <- match(sample_names, colnames(dfs[[i]]))
        dfs[[i]][, quant_cols_list[[i]]] <- dfs[[i]][, aligned_quant_cols, drop = FALSE]
        colnames(dfs[[i]])[quant_cols_list[[i]]] <- sample_names
        log_step("  Fraction '%s': reordered abundance columns to match '%s'.",
                 fraction_file_names[i], fraction_file_names[1])
      }
    }
  }

  for (i in seq_along(dfs)) {
    this_samples <- colnames(dfs[[i]])[quant_cols_list[[i]]]
    if (!identical(this_samples, sample_names)) {
      stop("Abundance columns must map to the same sample names/order across fractions after alignment. Mismatch in '",
           fraction_file_names[i], "'.")
    }
  }

  sample_info <- mp_build_sample_info(sample_meta_data, sample_names,
                                      "sample_meta_data", sample_id_column)
  resolved_sample_id_column <- attr(sample_info, "sample_col")

  mnar_variable_columns <- mp_resolve_mnar_columns(mnar_variables, sample_info,
                                                   resolved_sample_id_column)

  batch_groups <- NULL
  if (needs_batch_info) {
    batch_groups <- mp_build_batch_groups(sample_info, sample_names, batch_column)
  }
  if (!is.null(batch_specific_features_threshold) && length(batch_groups) < 2) {
    stop("`batch_specific_features_threshold` requires at least 2 distinct batch values in `batch_column`.")
  }
  if (batch_correction && length(batch_groups) < 2) {
    stop("`batch_column` must contain at least 2 distinct values for ComBat.")
  }

  log_step("Starting preprocessing: %d fraction(s), %d sample(s), %d total features.",
           length(dfs), length(sample_names), sum(vapply(dfs, nrow, integer(1))))

  if (zero_handling == "na") {
    log_step("Replacing zeros with NA.")
    for (i in seq_along(dfs)) {
      x <- dfs[[i]][, quant_cols_list[[i]], drop = FALSE]
      n_zero <- sum(x == 0, na.rm = TRUE)
      x[x == 0] <- NA
      dfs[[i]][, quant_cols_list[[i]]] <- x
      log_step("  Fraction '%s': replaced %d zero values.", fraction_file_names[i], n_zero)
    }
  } else {
    log_step("Keeping zero values unchanged (zero_handling='leave').")
  }

  cleaned <- lapply(seq_along(dfs), function(i) mp_remove_all_missing_rows(dfs[[i]], quant_cols_list[[i]]))
  for (i in seq_along(cleaned)) {
    dfs[[i]] <- cleaned[[i]]$data
    log_step("  Fraction '%s': removed %d all-missing features (kept %d/%d).",
             fraction_file_names[i], cleaned[[i]]$n_missing,
             cleaned[[i]]$n_after, cleaned[[i]]$n_before)
  }

  if (length(dfs) > 1) {
    for (i in seq_along(dfs)) {
      dfs[[i]][[precursor_id_column]] <- paste0(dfs[[i]][[precursor_id_column]], "_", fraction_file_names[i])
    }
  }

  se_list <- list()
  for (i in seq_along(dfs)) {
    se_list[[fraction_file_names[i]]] <- mp_run_quiet(
      QFeatures::readSummarizedExperiment(dfs[[i]],
                                          quantCols = quant_cols_list[[i]],
                                          fnames = precursor_id_column)
    )
    stopifnot(identical(
      colnames(SummarizedExperiment::assay(se_list[[fraction_file_names[i]]])),
      sample_names))
  }
  qf <- QFeatures::QFeatures(se_list)
  log_step("Created QFeatures object with assays: %s", paste(names(qf), collapse = ", "))

  # -------------------------------------------------------------------------
  # MNAR/MAR classification and filtering
  # -------------------------------------------------------------------------
  do_mnar <- (impute_method == "mixed")
  mnar_tables <- list()
  initial_assay_names <- names(qf)
  mar_thr  <- missing_feature_proportion_threshold[1]
  mnar_thr <- mp_mnar_threshold(missing_feature_proportion_threshold)

  if (do_mnar) {
    log_step("Starting MNAR classification for mixed imputation (p_adj <= %.4f).", mnar_significance_threshold)
    for (nm in initial_assay_names) {
      res <- mp_classify_missingness(qf, sample_info, mnar_variable_columns, nm,
                                     precursor_id_column, mnar_significance_threshold, mar_thr)
      qf <- res$object
      mnar_tables[[nm]] <- res$results
      log_step("  Assay '%s': MNAR features = %d / %d.", nm,
               sum(res$results$MNAR, na.rm = TRUE), nrow(res$results))
    }

    before_counts <- vapply(initial_assay_names, function(nm) nrow(qf[[nm]]), integer(1))
    qf <- mp_run_quiet(QFeatures::filterFeatures(
      qf,
      QFeatures::VariableFilter(field = "Remove_feature", value = 0, condition = "=="),
      i = initial_assay_names
    ))
    qf <- mp_run_quiet(QFeatures::filterFeatures(
      qf,
      QFeatures::VariableFilter(field = "Frac_missing", value = mnar_thr, condition = "<="),
      i = initial_assay_names
    ))
    after_counts <- vapply(initial_assay_names, function(nm) nrow(qf[[nm]]), integer(1))
    for (nm in initial_assay_names) {
      log_step("  Assay '%s': kept %d/%d features after MNAR-aware filtering.", nm, after_counts[nm], before_counts[nm])
    }
  } else {
    log_step("Applying MAR-only missingness filter at %.2f.", mar_thr)
    for (nm in names(qf)) {
      se_nm <- qf[[nm]]
      keep <- rowMeans(is.na(SummarizedExperiment::assay(se_nm))) <= mar_thr
      qf[[nm]] <- se_nm[keep, ]
      log_step("  Assay '%s': kept %d/%d features.", nm, sum(keep), length(keep))
    }
  }

  if (sum(vapply(names(qf), function(nm) nrow(qf[[nm]]), integer(1))) == 0) {
    stop("No features survived the missingness filter. Relax ",
         "`missing_feature_proportion_threshold` or check the input data.")
  }

  # -------------------------------------------------------------------------
  # Batch-specific missingness filter
  # -------------------------------------------------------------------------
  if (!is.null(batch_specific_features_threshold)) {
    log_step("Removing features with > %.2f missingness in any batch.", batch_specific_features_threshold)
    for (nm in names(qf)) {
      batch_filter <- mp_remove_batch_sparse_features(qf[[nm]], batch_groups, batch_specific_features_threshold)
      qf[[nm]] <- batch_filter$se
      log_step("  Assay '%s': removed %d batch-specific sparse features (kept %d/%d).",
               nm, batch_filter$n_removed, batch_filter$n_after, batch_filter$n_before)
      for (i in seq_len(nrow(batch_filter$details))) {
        log_step("    Batch '%s': %d features exceeded the threshold before union.",
                 batch_filter$details$Batch[i],
                 batch_filter$details$Features_above_threshold[i])
      }
    }
  } else {
    log_step("Skipping batch-specific missingness filter.")
  }

  # -------------------------------------------------------------------------
  # Imputation
  # -------------------------------------------------------------------------
  source_assays <- names(qf)
  if (impute_method != "none") {
    log_step("Starting imputation using method '%s'.", impute_method)
    n_imputation_steps <- if (use_batchwise_imputation) {
      length(source_assays) * length(batch_groups)
    } else {
      length(source_assays)
    }
    resolved_imputation_seeds <- mp_resolve_seeds(imputation_seeds, n_imputation_steps)
    next_seed_idx <- 1L
    for (nm in source_assays) {
      na_before <- sum(is.na(SummarizedExperiment::assay(qf[[nm]])))
      if (use_batchwise_imputation) {
        if (impute_method == "mixed") {
          batch_all_missing <- mp_find_batch_all_missing_features(qf[[nm]], batch_groups)
          if (nrow(batch_all_missing) > 0) {
            batch_counts <- sort(table(batch_all_missing$Batch), decreasing = TRUE)
            stop(
              "Batch-wise mixed imputation cannot be performed because ",
              nrow(batch_all_missing), " feature(s) are completely missing within at least one batch (",
              paste(sprintf("%s=%d", names(batch_counts), as.integer(batch_counts)), collapse = ", "),
              "). Example feature IDs: ",
              paste(utils::head(unique(batch_all_missing$Feature), 5), collapse = ", "),
              ". Use a stricter `missing_feature_proportion_threshold`, provide `batch_specific_features_threshold`, or set `impute_behaviour = 'global'`."
            )
          }
        }
        log_step("  Assay '%s': imputing missing values separately within %d batch(es).", nm, length(batch_groups))
        batch_seeds <- resolved_imputation_seeds[next_seed_idx:(next_seed_idx + length(batch_groups) - 1L)]
        next_seed_idx <- next_seed_idx + length(batch_groups)
        imp <- mp_impute_assay_by_batch(qf[[nm]], batch_groups, impute_method,
                                        seeds = batch_seeds,
                                        impute_with_value = impute_with_value,
                                        verbose = verbose)
      } else {
        step_seed <- resolved_imputation_seeds[next_seed_idx]
        next_seed_idx <- next_seed_idx + 1L
        imp <- mp_impute_single_assay(qf[[nm]], impute_method, paste0("Assay '", nm, "'"),
                                      seed = step_seed,
                                      impute_with_value = impute_with_value,
                                      verbose = verbose)
      }

      na_after <- sum(is.na(SummarizedExperiment::assay(imp)))
      if (na_after > 0) {
        stop("Imputation method '", impute_method, "' returned ", na_after,
             " missing values in assay '", nm, "'.")
      }
      mp_assert_finite(SummarizedExperiment::assay(imp), "after imputation")

      qf[[paste0(nm, "_imputed")]] <- imp
      log_step("  Assay '%s': imputation complete (%d -> %d missing values).", nm, na_before, na_after)
    }
    use_assays <- paste0(source_assays, "_imputed")
  } else {
    log_step("Skipping imputation (`impute_method='none'`).")
    use_assays <- source_assays
  }

  # -------------------------------------------------------------------------
  # Combine fractions
  # -------------------------------------------------------------------------
  if (length(use_assays) > 1) {
    log_step("Combining %d assays into a single raw assay.", length(use_assays))
    # bind_rows takes the union of annotation columns; base rbind() would
    # reject fractions that carry different ancillary annotations, which the
    # documented contract permits.
    combined <- dplyr::bind_rows(lapply(use_assays, function(nm) {
      cbind(as.data.frame(SummarizedExperiment::rowData(qf[[nm]])),
            SummarizedExperiment::assay(qf, nm))
    }))
    mp_validate_feature_ids(combined[[precursor_id_column]], precursor_id_column,
                            "the combined fractions")
    quant_idx <- match(sample_names, colnames(combined))
    stopifnot(!anyNA(quant_idx))
    se_combined <- mp_run_quiet(QFeatures::readSummarizedExperiment(
      combined,
      quantCols = quant_idx,
      fnames = precursor_id_column
    ))
    stopifnot(identical(colnames(SummarizedExperiment::assay(se_combined)), sample_names))
    q <- QFeatures::QFeatures(list(raw = se_combined))
  } else {
    q <- QFeatures::QFeatures(list(raw = qf[[use_assays[1]]]))
  }
  assay_name <- "raw"

  # -------------------------------------------------------------------------
  # QC plots
  # -------------------------------------------------------------------------
  plots <- list()
  if (plot_qc_charts) {
    plots$raw_quantities <- mp_plot_density(SummarizedExperiment::assay(q[[assay_name]]), "Raw", show_plot_legends)
  }

  # -------------------------------------------------------------------------
  # Transform + normalise
  # -------------------------------------------------------------------------
  normalised <- mp_transform_and_normalise(q, assay_name, normalisation_method,
                                           plot_qc_charts, show_plot_legends, verbose)
  q <- normalised$q
  assay_name <- normalised$assay_name
  plots <- c(plots, normalised$plots)

  # -------------------------------------------------------------------------
  # Batch correction
  # -------------------------------------------------------------------------
  if (batch_correction) {
    log_step("Applying ComBat batch correction with batch column '%s'.", batch_column)
    x <- SummarizedExperiment::assay(q[[assay_name]])

    if (plot_qc_charts) {
      before_pca_plot <- mp_build_batch_pca_plot(
        assay_matrix = x, sample_info = sample_info, batch_col = batch_column,
        title = "Before batch correction", context_label = "pre-batch-correction",
        verbose = verbose)
      if (!is.null(before_pca_plot)) plots$before_batch_correction <- before_pca_plot
    }

    corrected <- mp_run_combat(x, sample_info, batch_column,
                               batch_correction_variable_columns, verbose)
    batch_corrected_se <- q[[assay_name]]
    SummarizedExperiment::assay(batch_corrected_se) <- corrected
    q[["batch_corrected"]] <- batch_corrected_se
    assay_name <- "batch_corrected"
    log_step("Batch correction complete.")

    if (plot_qc_charts) {
      plots$batch_corrected_quantities <- mp_plot_density(
        SummarizedExperiment::assay(q[[assay_name]]), "After batch correction", show_plot_legends)
      after_pca_plot <- mp_build_batch_pca_plot(
        assay_matrix = SummarizedExperiment::assay(q[[assay_name]]),
        sample_info = sample_info, batch_col = batch_column,
        title = "After batch correction", context_label = "post-batch-correction",
        verbose = verbose)
      if (!is.null(after_pca_plot)) plots$after_batch_correction <- after_pca_plot
    }
  } else {
    log_step("Batch correction skipped.")
  }

  # -------------------------------------------------------------------------
  # Aggregate to proteins
  # -------------------------------------------------------------------------
  aggregation_performed <- !is.null(protein_aggregator_column)
  final_assay_name <- assay_name
  abundances <- as.data.frame(SummarizedExperiment::assay(q[[final_assay_name]]))

  if (aggregation_performed) {
    row_data <- SummarizedExperiment::rowData(q[[assay_name]])
    if (!(protein_aggregator_column %in% colnames(row_data))) {
      stop("`protein_aggregator_column` ('", protein_aggregator_column,
           "') not found in rowData of assay '", assay_name, "'.")
    }

    protein_ids <- as.character(row_data[[protein_aggregator_column]])
    if (anyNA(protein_ids) || any(!nzchar(protein_ids))) {
      stop("`protein_aggregator_column` contains missing or empty protein identifiers after preprocessing.")
    }

    log_step("Aggregating precursor features to proteins using '%s'.", protein_aggregator_column)
    q <- mp_run_quiet(QFeatures::aggregateFeatures(
      q,
      i = assay_name,
      fcol = protein_aggregator_column,
      name = "aggregated",
      fun = protein_aggregator_method,
      na.rm = TRUE
    ))

    proteins_mat <- SummarizedExperiment::assay(q[["aggregated"]])
    proteins_na <- sum(is.na(proteins_mat))
    if (proteins_na > 0) {
      stop("Protein aggregation produced ", proteins_na, " missing values.")
    }
    mp_assert_finite(proteins_mat, "after aggregation")
    if (plot_qc_charts) {
      plots$aggregated_quantities <- mp_plot_density(proteins_mat, "After aggregation", show_plot_legends)
    }
    final_assay_name <- "aggregated"
    abundances <- as.data.frame(proteins_mat)
    log_step("Finished: %d proteins x %d samples.", nrow(abundances), ncol(abundances))
  } else {
    log_step("Aggregation skipped (`protein_aggregator_column = NULL`).")
    log_step("Finished: %d features x %d samples.", nrow(abundances), ncol(abundances))
  }

  out <- list(
    abundances = abundances,
    final_assay_name = final_assay_name,
    qfeatures = q
  )
  if (aggregation_performed) {
    out$proteins <- abundances
  } else {
    out$precursors <- abundances
  }
  if (return_mnar_results) {
    if (do_mnar && length(mnar_tables) == 0) {
      stop("`mnar_results` is empty while `impute_method = 'mixed'`. MNAR classification failed.")
    }
    out$mnar_results <- mnar_tables
  }
  if (plot_qc_charts) out$plots <- plots
  out
}
