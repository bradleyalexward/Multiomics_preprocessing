# Internal helpers shared by prepare_metabolites() and prepare_proteins().
#
# Nothing in this file is exported. Both public functions previously carried
# their own copy of most of these; keeping one definition means a fix cannot
# be applied to one modality and forgotten in the other.

utils::globalVariables(c("PC1", "PC2", "Batch"))

# --- small utilities -------------------------------------------------------

mp_log <- function(verbose, fmt, ...) {
  if (isTRUE(verbose)) message(sprintf(fmt, ...))
  invisible(NULL)
}

# Runs `expr` with printed output and progress messages hidden. Warnings are
# only hidden when explicitly requested: vsn convergence problems and ComBat
# diagnostics are worth seeing.
mp_run_quiet <- function(expr, suppress_warnings = TRUE) {
  value <- NULL
  if (suppress_warnings) {
    suppressWarnings(suppressMessages(
      utils::capture.output(value <- expr, type = "output")))
  } else {
    suppressMessages(
      utils::capture.output(value <- expr, type = "output"))
  }
  value
}

mp_assert_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

mp_assert_string <- function(x, arg, allow_null = FALSE) {
  if (allow_null && is.null(x)) return(invisible(TRUE))
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop("`", arg, "` must be a single non-empty character string.", call. = FALSE)
  }
  invisible(TRUE)
}

# --- argument validation shared by both workflows ---------------------------

mp_validate_thresholds <- function(missing_feature_proportion_threshold,
                                   batch_specific_features_threshold,
                                   batch_column_missing) {
  x <- missing_feature_proportion_threshold
  if (!is.numeric(x) || !(length(x) %in% c(1, 2)) || anyNA(x) ||
      any(x < 0 | x > 1)) {
    stop("`missing_feature_proportion_threshold` must be numeric length 1 or 2 ",
         "with values in [0, 1].", call. = FALSE)
  }
  if (length(x) == 2 && x[2] < x[1]) {
    stop("When length 2, `missing_feature_proportion_threshold` must be ",
         "c(MAR_threshold, MNAR_threshold) with MNAR >= MAR.", call. = FALSE)
  }

  if (!is.null(batch_specific_features_threshold)) {
    b <- batch_specific_features_threshold
    if (!is.numeric(b) || length(b) != 1 || is.na(b) || b < 0 || b > 1) {
      stop("`batch_specific_features_threshold` must be `NULL` or a single ",
           "numeric value in [0, 1].", call. = FALSE)
    }
    if (isTRUE(batch_column_missing)) {
      stop("When `batch_specific_features_threshold` is supplied, ",
           "`batch_column` must also be provided.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Returns the resolved MNAR threshold. A single value sets the MAR threshold;
# MNAR features keep the documented 0.75 allowance, but never a stricter one
# than MAR features get, which would invert the policy the two-value form
# explicitly forbids.
mp_mnar_threshold <- function(missing_feature_proportion_threshold) {
  if (length(missing_feature_proportion_threshold) >= 2) {
    missing_feature_proportion_threshold[2]
  } else {
    max(0.75, missing_feature_proportion_threshold[1])
  }
}

mp_validate_normalisation <- function(normalisation_method) {
  if (!is.character(normalisation_method)) {
    stop("`normalisation_method` must be a character vector.", call. = FALSE)
  }
  known <- MsCoreUtils::normalizeMethods()
  if (length(normalisation_method) == 1) {
    if (!(normalisation_method[1] %in% known)) {
      stop("Unsupported normalisation method: '", normalisation_method[1], "'.",
           call. = FALSE)
    }
    mp_check_normalisation_engine(normalisation_method[1])
  } else if (length(normalisation_method) == 2) {
    if (!(normalisation_method[1] %in% c("log2", "log", "log10"))) {
      stop("When length 2, first normalisation element must be one of: ",
           "log2, log, log10.", call. = FALSE)
    }
    if (!(normalisation_method[2] %in% known)) {
      stop("Unsupported normalisation method: '", normalisation_method[2], "'.",
           call. = FALSE)
    }
    mp_check_normalisation_engine(normalisation_method[2])
  } else {
    stop("`normalisation_method` must have length 1 or 2.", call. = FALSE)
  }
  invisible(TRUE)
}

mp_validate_mnar_variables <- function(mnar_variables, impute_method) {
  if (impute_method == "mixed") {
    if (is.null(mnar_variables) || !is.character(mnar_variables) ||
        length(mnar_variables) == 0 ||
        any(is.na(mnar_variables) | mnar_variables == "")) {
      stop("`mnar_variables` must be supplied as one or more column names from ",
           "`sample_meta_data` when `impute_method = 'mixed'`.", call. = FALSE)
    }
  } else if (!is.null(mnar_variables) &&
             (!is.character(mnar_variables) ||
              any(is.na(mnar_variables) | mnar_variables == ""))) {
    stop("`mnar_variables` must be `NULL` or a character vector of column names ",
         "from `sample_meta_data`.", call. = FALSE)
  }
  invisible(TRUE)
}

# Resolves mnar_variables against the metadata and rejects the sample ID column.
mp_resolve_mnar_columns <- function(mnar_variables, sample_info, sample_id_column) {
  if (is.null(mnar_variables)) return(NULL)
  cols <- unique(as.character(mnar_variables))
  missing_cols <- setdiff(cols, colnames(sample_info))
  if (length(missing_cols) > 0) {
    stop("`mnar_variables` contains columns not found in `sample_meta_data`: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  forbidden <- intersect(cols, unique(c("Sample", sample_id_column)))
  if (length(forbidden) > 0) {
    stop("`mnar_variables` must reference covariate columns in ",
         "`sample_meta_data`, not the sample ID column: ",
         paste(forbidden, collapse = ", "), call. = FALSE)
  }
  cols
}

mp_validate_dataset_names <- function(names, n_files, arg) {
  if (is.null(names)) return(invisible(TRUE))
  if (!is.character(names)) {
    stop("`", arg, "` must be a character vector when supplied.", call. = FALSE)
  }
  if (length(names) != n_files) {
    stop("`", arg, "` must have one element per input file.", call. = FALSE)
  }
  if (any(is.na(names) | names == "")) {
    stop("`", arg, "` cannot contain NA or empty values.", call. = FALSE)
  }
  if (any(duplicated(names))) {
    stop("`", arg, "` must be unique.", call. = FALSE)
  }
  invisible(TRUE)
}

mp_validate_common_args <- function(input_files, sample_meta_data, verbose,
                                    batch_correction, plot_qc_charts,
                                    return_mnar_results,
                                    mnar_significance_threshold,
                                    impute_method, impute_with_value) {
  mp_assert_flag(verbose, "verbose")
  mp_assert_flag(batch_correction, "batch_correction")
  mp_assert_flag(plot_qc_charts, "plot_qc_charts")
  mp_assert_flag(return_mnar_results, "return_mnar_results")

  if (!is.list(input_files) || inherits(input_files, "data.frame") ||
      length(input_files) == 0) {
    stop("`input_files` must be a non-empty list of one or more data.frames.",
         call. = FALSE)
  }
  if (!all(vapply(input_files, inherits, logical(1), what = "data.frame"))) {
    stop("All elements of `input_files` must be data.frames.", call. = FALSE)
  }
  if (!is.data.frame(sample_meta_data) || nrow(sample_meta_data) == 0) {
    stop("`sample_meta_data` must be a non-empty data.frame.", call. = FALSE)
  }
  if (!is.numeric(mnar_significance_threshold) ||
      length(mnar_significance_threshold) != 1 ||
      is.na(mnar_significance_threshold) ||
      mnar_significance_threshold <= 0 || mnar_significance_threshold > 1) {
    stop("`mnar_significance_threshold` must be a single numeric value in (0, 1].",
         call. = FALSE)
  }
  if (impute_method == "MLE2") {
    stop("`impute_method = 'MLE2'` is defunct in QFeatures. Please use ",
         "`impute_method = 'MLE'`.", call. = FALSE)
  }
  if (impute_method == "with") {
    if (is.null(impute_with_value) || !is.numeric(impute_with_value) ||
        length(impute_with_value) != 1 || !is.finite(impute_with_value)) {
      stop("When `impute_method = 'with'`, `impute_with_value` must be a single ",
           "finite numeric value.", call. = FALSE)
    }
  }
  mp_check_impute_engines(impute_method)
  invisible(TRUE)
}

mp_validate_batch_args <- function(needs_batch_info, batch_column,
                                   batch_correction,
                                   batch_correction_variable_columns,
                                   sample_meta_data) {
  if (needs_batch_info) {
    mp_assert_string(batch_column, "batch_column")
    if (!(batch_column %in% colnames(sample_meta_data))) {
      stop("`batch_column` ('", batch_column, "') is not present in ",
           "`sample_meta_data`.", call. = FALSE)
    }
  }
  if (batch_correction && !is.null(batch_correction_variable_columns)) {
    if (!is.character(batch_correction_variable_columns)) {
      stop("`batch_correction_variable_columns` must be a character vector when ",
           "supplied.", call. = FALSE)
    }
    missing_cols <- setdiff(batch_correction_variable_columns,
                            colnames(sample_meta_data))
    if (length(missing_cols) > 0) {
      stop("Missing `batch_correction_variable_columns` in `sample_meta_data`: ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
  }
  invisible(TRUE)
}

mp_validate_aggregator <- function(column, method, arg_column, arg_method) {
  if (is.null(column)) return(invisible(TRUE))
  if (!is.function(method)) {
    stop("`", arg_method, "` must be a function when aggregation is requested.",
         call. = FALSE)
  }
  mp_assert_string(column, arg_column)
  invisible(TRUE)
}

# --- random number generator protection ------------------------------------

# The imputation steps switch to the L'Ecuyer-CMRG generator, which is the
# right choice for reproducible parallel work but must not leak into the
# caller's session. Save on entry, restore on exit (including on error).
mp_save_rng <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

mp_restore_rng <- function(old_seed) {
  if (is.null(old_seed)) {
    suppressWarnings(rm(".Random.seed", envir = globalenv()))
  } else {
    assign(".Random.seed", old_seed, envir = globalenv())
  }
  invisible(NULL)
}

# --- seeds -----------------------------------------------------------------

mp_validate_seeds <- function(seeds) {
  if (is.null(seeds)) return(NULL)
  if (!is.numeric(seeds) || length(seeds) == 0 || anyNA(seeds) ||
      any(!is.finite(seeds)) || any(seeds != floor(seeds))) {
    stop("`imputation_seeds` must be NULL or an integer vector.", call. = FALSE)
  }
  as.integer(seeds)
}

mp_resolve_seeds <- function(seeds, n_steps) {
  if (n_steps < 1L) return(integer(0))
  if (is.null(seeds)) return(rep(NA_integer_, n_steps))
  if (length(seeds) == 1L) return(as.integer(seeds[1] + seq_len(n_steps) - 1L))
  if (length(seeds) == n_steps) return(seeds)
  stop("`imputation_seeds` must have length 1 or ", n_steps,
       " for the current imputation configuration.", call. = FALSE)
}

# --- optional imputation engines -------------------------------------------

# MsCoreUtils implements several methods by calling packages it only Suggests,
# so they are not installed with this package. Fail here, naming the package,
# rather than several minutes later from inside Bioconductor.
mp_impute_engines <- function(method) {
  switch(method,
         mixed   = c("missForest", "imputeLCMD"),
         RF      = "missForest",
         QRILC   = ,
         MinDet  = ,
         MinProb = "imputeLCMD",
         knn     = "impute",
         bpca    = "pcaMethods",
         MLE     = "norm",
         character(0))
}

mp_check_impute_engines <- function(method) {
  needed <- mp_impute_engines(method)
  if (length(needed) == 0) return(invisible(TRUE))
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent) > 0) {
    stop("`impute_method = '", method, "'` requires the ",
         paste(absent, collapse = " and "), " package",
         if (length(absent) > 1) "s" else "", ", which ",
         if (length(absent) > 1) "are" else "is", " not installed. Install with:\n",
         "  BiocManager::install(c(\"", paste(absent, collapse = "\", \""), "\"))",
         call. = FALSE)
  }
  invisible(TRUE)
}

mp_check_normalisation_engine <- function(method) {
  if (method %in% c("quantiles", "quantiles.robust") &&
      !requireNamespace("preprocessCore", quietly = TRUE)) {
    stop("`normalisation_method = '", method, "'` requires the preprocessCore ",
         "package, which is not installed. Install with:\n",
         "  BiocManager::install(\"preprocessCore\")", call. = FALSE)
  }
  if (method == "vsn") {
    # vsn is called by MsCoreUtils, which only suggests it. It is the default
    # normalisation here, so it is a hard dependency of this package instead.
    # This reference is what declares that to R CMD check -- it is load-bearing,
    # not a leftover.
    invisible(vsn::justvsn)
  }
  invisible(TRUE)
}

# --- numerical safety gates ------------------------------------------------

mp_assert_finite <- function(x, stage) {
  bad <- !is.finite(x) & !is.na(x)
  if (any(bad)) {
    stop("Non-finite values (NaN or Inf) were produced ", stage, " in ",
         sum(bad), " cell(s). This usually follows a transform applied to ",
         "zero or negative abundances.", call. = FALSE)
  }
  invisible(TRUE)
}

# log(0) is -Inf and log(negative) is undefined. Whether a zero means "truly
# absent" or "not measured" is a scientific judgement, so ask rather than
# silently adding a pseudocount.
mp_assert_log_domain <- function(x, transform) {
  bad <- !is.na(x) & x <= 0
  if (any(bad)) {
    stop("`normalisation_method` requests a ", transform, " transform, but ",
         sum(bad), " abundance value(s) are zero or negative. Set ",
         "`zero_handling = 'na'` and impute, choose a non-logarithmic ",
         "normalisation, or transform the data upstream.", call. = FALSE)
  }
  invisible(TRUE)
}

# --- column and sample resolution ------------------------------------------

# Numeric indices and character names both preserve the order the user asked
# for. A logical mask cannot, so it is converted to positions here and the
# assay is always built from positions.
mp_resolve_quant_cols <- function(df, cols, label, arg_name) {
  if (is.numeric(cols)) {
    idx <- unique(as.integer(cols))
    if (anyNA(idx) || any(idx < 1 | idx > ncol(df))) {
      stop("`", arg_name, "` contains invalid numeric indices for '", label, "'.",
           call. = FALSE)
    }
    return(idx)
  }
  if (is.character(cols)) {
    if (anyDuplicated(cols) > 0) {
      stop("`", arg_name, "` contains duplicated column names: ",
           paste(unique(cols[duplicated(cols)]), collapse = ", "), ".", call. = FALSE)
    }
    missing_cols <- setdiff(cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop("`", arg_name, "` contains missing column names for '", label, "': ",
           paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    return(match(cols, colnames(df)))
  }
  if (is.logical(cols)) {
    if (length(cols) != ncol(df)) {
      stop("Logical `", arg_name, "` must have length ncol(df) for '", label, "'.",
           call. = FALSE)
    }
    idx <- which(cols)
    if (length(idx) == 0) {
      stop("Logical `", arg_name, "` selected no columns for '", label, "'.",
           call. = FALSE)
    }
    return(idx)
  }
  stop("`", arg_name, "` must be numeric, character, or logical.", call. = FALSE)
}

# as.numeric() on a factor returns the internal level codes, so a column
# displaying 10, 100, 20 would silently become 1, 2, 3. Go via the labels.
mp_coerce_abundance_to_numeric <- function(df, quant_cols, label) {
  for (j in quant_cols) {
    value <- df[[j]]
    if (is.numeric(value)) {
      if (any(!is.na(value) & !is.finite(value))) {
        stop("Non-finite values found in abundance column '", colnames(df)[j],
             "' for '", label, "'.", call. = FALSE)
      }
      next
    }
    if (is.factor(value)) value <- as.character(value)
    converted <- suppressWarnings(as.numeric(value))
    if (any(!is.na(value) & is.na(converted))) {
      stop("Non-numeric values found in abundance column '", colnames(df)[j],
           "' for '", label, "'.", call. = FALSE)
    }
    if (any(!is.na(converted) & !is.finite(converted))) {
      stop("Non-finite values found in abundance column '", colnames(df)[j],
           "' for '", label, "'.", call. = FALSE)
    }
    df[[j]] <- converted
  }
  df
}

# which.max() silently returns the first winner, so two metadata columns that
# both hold every sample ID would resolve to whichever came first. Refuse.
mp_detect_sample_col <- function(df, assay_sample_names, label) {
  matches <- vapply(df, function(x) sum(as.character(x) %in% assay_sample_names),
                    numeric(1))
  full_cover <- names(matches)[matches == length(assay_sample_names)]
  if (length(full_cover) == 0) {
    stop("Could not identify a sample ID column in ", label,
         " by matching assay sample names. Set `sample_id_column` explicitly.",
         call. = FALSE)
  }
  if (length(full_cover) > 1) {
    stop("Multiple columns in ", label, " contain all assay sample IDs: ",
         paste(full_cover, collapse = ", "),
         ". Set `sample_id_column` to choose one.", call. = FALSE)
  }
  full_cover
}

# Always returns a base data.frame with an explicit Sample column. Tibbles do
# not carry row names through dplyr verbs, which used to silently break the
# QC plots; an ordinary column is a key the data structure actually keeps.
mp_build_sample_info <- function(df, assay_sample_names, label,
                                 sample_id_column = NULL) {
  meta <- as.data.frame(df, stringsAsFactors = FALSE)

  if (is.null(sample_id_column)) {
    sample_col <- mp_detect_sample_col(meta, assay_sample_names, label)
  } else {
    mp_assert_string(sample_id_column, "sample_id_column")
    if (!(sample_id_column %in% colnames(meta))) {
      stop("`sample_id_column` ('", sample_id_column, "') is not a column of ",
           label, ".", call. = FALSE)
    }
    sample_col <- sample_id_column
  }

  meta$Sample <- as.character(meta[[sample_col]])
  if (anyNA(meta$Sample) || any(!nzchar(meta$Sample))) {
    stop(label, " contains missing or empty sample identifiers in column '",
         sample_col, "'.", call. = FALSE)
  }
  if (anyDuplicated(meta$Sample) > 0) {
    stop(label, " contains duplicated sample identifiers in column '",
         sample_col, "'.", call. = FALSE)
  }
  missing_samples <- setdiff(assay_sample_names, meta$Sample)
  if (length(missing_samples) > 0) {
    stop(label, " is missing sample IDs required by the assay: ",
         paste(missing_samples, collapse = ", "), call. = FALSE)
  }

  meta <- meta[match(assay_sample_names, meta$Sample), , drop = FALSE]
  rownames(meta) <- meta$Sample
  attr(meta, "sample_col") <- sample_col
  meta
}

# Reorders metadata to the assay's own column order and proves it matched.
# Every positional handover (ComBat, the batch vector, the design matrix)
# goes through here, so a misalignment stops rather than silently relabelling.
mp_align_metadata <- function(sample_info, assay_colnames, cols) {
  idx <- match(assay_colnames, sample_info$Sample)
  if (anyNA(idx)) {
    stop("Assay columns are absent from the sample metadata: ",
         paste(assay_colnames[is.na(idx)], collapse = ", "), call. = FALSE)
  }
  aligned <- sample_info[idx, unique(c("Sample", cols)), drop = FALSE]
  stopifnot(identical(as.character(aligned$Sample), as.character(assay_colnames)))
  aligned
}

mp_validate_feature_ids <- function(ids, id_col, label) {
  ids <- as.character(ids)
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop("`", id_col, "` contains missing, empty, or whitespace-only values in '",
         label, "'.", call. = FALSE)
  }
  dups <- unique(ids[duplicated(ids)])
  if (length(dups) > 0) {
    stop("`", id_col, "` must uniquely identify rows within '", label,
         "'. Duplicated: ", paste(utils::head(dups, 5), collapse = ", "),
         if (length(dups) > 5) sprintf(" (and %d more)", length(dups) - 5) else "",
         ".", call. = FALSE)
  }
  invisible(TRUE)
}

mp_remove_all_missing_rows <- function(df, quant_cols) {
  x <- df[, quant_cols, drop = FALSE]
  all_missing <- rowSums(is.na(x)) == ncol(x)
  list(data      = df[!all_missing, , drop = FALSE],
       n_missing = sum(all_missing),
       n_before  = nrow(df),
       n_after   = sum(!all_missing))
}

# --- missingness classification --------------------------------------------

# For each feature, test whether detection (observed vs not) is explained by
# the chosen sample covariates.
#
# This replaces an earlier implementation that took the smallest Wald p-value
# across all model coefficients *including the intercept*. The intercept only
# asks whether a feature's overall detection rate differs from 50%, which is
# not the question, and the Wald test collapses towards p = 1 under complete
# separation -- exactly the case where the covariate explains detection
# perfectly, and therefore the clearest signal there is.
#
# A single likelihood-ratio test comparing the covariate model against an
# intercept-only model asks the intended question in one number and is not
# affected by separation. glm.fit() reports both deviances from one fit, so
# no second model is needed.
#
# Note this detects association with *observed* metadata. That is a practical
# heuristic for routing features to left-censored imputation, not a proof of
# the formal missing-not-at-random mechanism.
mp_classify_missingness <- function(object, sample_info, covariate_cols,
                                    assay_name, feature_name, p_adj_cutoff,
                                    mar_threshold) {
  assay_mat <- SummarizedExperiment::assay(object[[assay_name]])
  detected  <- !is.na(assay_mat)

  # Covariates are sample-level, so the design matrix is identical for every
  # feature. Build it once instead of rebuilding it inside the loop.
  covars <- mp_align_metadata(sample_info, colnames(detected), covariate_cols)
  if (anyNA(covars[, covariate_cols, drop = FALSE])) {
    stop("`mnar_variables` contains missing values for one or more assay ",
         "samples. Complete those columns or drop the affected samples.",
         call. = FALSE)
  }
  mm <- stats::model.matrix(stats::reformulate(covariate_cols), data = covars)
  if (ncol(mm) < 2L) {
    stop("`mnar_variables` (", paste(covariate_cols, collapse = ", "),
         ") does not vary across the assay samples, so there is nothing to ",
         "test. Choose a covariate with at least two levels.", call. = FALSE)
  }

  n_feat    <- nrow(detected)
  pvals     <- rep(1, n_feat)
  converged <- rep(NA, n_feat)

  for (i in seq_len(n_feat)) {
    fit <- suppressWarnings(try(
      stats::glm.fit(mm, as.numeric(detected[i, ]), family = stats::binomial()),
      silent = TRUE))
    if (inherits(fit, "try-error")) next
    stat <- fit$null.deviance - fit$deviance
    df   <- fit$df.null - fit$df.residual
    if (is.finite(stat) && df > 0) {
      pvals[i] <- stats::pchisq(max(stat, 0), df, lower.tail = FALSE)
    }
    converged[i] <- isTRUE(fit$converged)
  }

  # One multiple-testing correction, across features. The previous code also
  # ran a correction within each feature before taking the minimum.
  adj_pvals    <- stats::p.adjust(pvals, method = "BH")
  frac_missing <- rowMeans(!detected)

  results <- data.frame(
    Feature         = rownames(detected),
    pvals           = pvals,
    adj.pvals       = adj_pvals,
    MNAR            = adj_pvals <= p_adj_cutoff,
    Frac_missing    = frac_missing,
    Model_converged = converged,
    Remove_feature  = (adj_pvals > p_adj_cutoff) & (frac_missing > mar_threshold),
    row.names        = NULL,
    stringsAsFactors = FALSE
  )
  colnames(results)[1] <- feature_name

  rd <- SummarizedExperiment::rowData(object[[assay_name]])
  idx <- match(rownames(object[[assay_name]]), results[[feature_name]])
  stopifnot(!anyNA(idx))
  for (col in setdiff(colnames(results), feature_name)) {
    rd[[col]] <- results[[col]][idx]
  }
  SummarizedExperiment::rowData(object[[assay_name]]) <- rd

  list(object = object, results = results)
}

# --- batches ---------------------------------------------------------------

mp_build_batch_groups <- function(sample_info, assay_sample_names, batch_col) {
  if (!(batch_col %in% colnames(sample_info))) {
    stop("`batch_column` ('", batch_col, "') is not present in `sample_meta_data`.",
         call. = FALSE)
  }
  aligned <- mp_align_metadata(sample_info, assay_sample_names, batch_col)
  batch_values <- as.character(aligned[[batch_col]])
  if (anyNA(batch_values) || any(!nzchar(batch_values))) {
    stop("`batch_column` contains missing or empty values for one or more ",
         "assay samples.", call. = FALSE)
  }
  split(assay_sample_names, batch_values)
}

mp_remove_batch_sparse_features <- function(se_obj, batch_groups, threshold) {
  x <- SummarizedExperiment::assay(se_obj)
  per_batch <- vapply(batch_groups,
                      function(s) rowMeans(is.na(x[, s, drop = FALSE])) > threshold,
                      logical(nrow(x)))
  per_batch <- matrix(per_batch, nrow = nrow(x),
                      dimnames = list(NULL, names(batch_groups)))
  remove_mask <- rowSums(per_batch) > 0
  list(se       = se_obj[!remove_mask, ],
       details  = data.frame(Batch = names(batch_groups),
                             Features_above_threshold = colSums(per_batch),
                             row.names = NULL, stringsAsFactors = FALSE),
       n_removed = sum(remove_mask),
       n_before  = nrow(x),
       n_after   = sum(!remove_mask))
}

mp_find_batch_all_missing_features <- function(se_obj, batch_groups) {
  x <- SummarizedExperiment::assay(se_obj)
  ids <- rownames(x)
  if (is.null(ids)) ids <- as.character(seq_len(nrow(x)))

  hits <- lapply(names(batch_groups), function(b) {
    all_missing <- rowSums(!is.na(x[, batch_groups[[b]], drop = FALSE])) == 0
    if (!any(all_missing)) return(NULL)
    data.frame(Feature = ids[all_missing], Batch = b,
               row.names = NULL, stringsAsFactors = FALSE)
  })
  hits <- Filter(Negate(is.null), hits)
  if (length(hits) == 0) {
    return(data.frame(Feature = character(0), Batch = character(0),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, hits)
}

# --- imputation ------------------------------------------------------------

# Mixed imputation splits the matrix in two and imputes each half separately:
# QRILC for the candidate-MNAR features, random forest for the rest. Both
# estimate a distribution per sample, so a block in which some sample has fewer
# than two observed values cannot be imputed. Left unchecked this surfaces as
# "'upper' not specified or contains NA" from tmvtnorm, three packages away
# from anything the user wrote.
mp_check_mixed_blocks <- function(x, randna, label, min_observed = 10L) {
  blocks <- list(`random forest` = x[randna, , drop = FALSE],
                 QRILC           = x[!randna, , drop = FALSE])
  for (nm in names(blocks)) {
    block <- blocks[[nm]]
    if (nrow(block) == 0) next
    observed <- colSums(!is.na(block))

    # Fewer than two observed values makes the per-sample fit undefined.
    impossible <- names(observed)[observed < 2L]
    if (length(impossible) > 0) {
      stop("Mixed imputation cannot proceed for ", label, ": the ", nm,
           " block holds ", nrow(block), " feature(s), leaving ",
           length(impossible), " sample(s) with fewer than two observed values (",
           paste(utils::head(impossible, 5), collapse = ", "),
           if (length(impossible) > 5) ", ..." else "",
           "). There is no distribution to estimate. Either choose a single ",
           "`impute_method` instead of 'mixed', relax ",
           "`mnar_significance_threshold` so more features are classified as ",
           "MNAR, or use `impute_behaviour = 'global'` if splitting by batch ",
           "made the block too small.", call. = FALSE)
    }

    # Above that it will run, but a handful of points is a poor basis for a
    # per-sample distribution, so say so rather than returning quiet nonsense.
    thin <- names(observed)[observed < min_observed]
    if (length(thin) > 0) {
      warning("Mixed imputation for ", label, ": the ", nm, " block gives ",
              length(thin), " sample(s) fewer than ", min_observed,
              " observed values (", paste(utils::head(thin, 5), collapse = ", "),
              if (length(thin) > 5) ", ..." else "",
              "). Imputed values for those samples rest on a thin estimate.",
              call. = FALSE)
    }
  }
  invisible(TRUE)
}

mp_impute_single_assay <- function(se_obj, method, label, seed = NA_integer_,
                                   impute_with_value = NULL, verbose = TRUE) {
  if (!is.na(seed)) {
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    mp_log(verbose, "  %s: using imputation seed %d.", label, seed)
  }
  if (method == "mixed") {
    randna <- !as.logical(SummarizedExperiment::rowData(se_obj)$MNAR)
    randna[is.na(randna)] <- TRUE
    mp_check_mixed_blocks(SummarizedExperiment::assay(se_obj), randna, label)
    mp_run_quiet(QFeatures::impute(se_obj, method = "mixed", randna = randna,
                                   mar = "RF", mnar = "QRILC"))
  } else if (method == "with") {
    mp_run_quiet(QFeatures::impute(se_obj, method = "with", val = impute_with_value))
  } else {
    mp_run_quiet(QFeatures::impute(se_obj, method = method))
  }
}

mp_impute_assay_by_batch <- function(se_obj, batch_groups, method, seeds,
                                     impute_with_value = NULL, verbose = TRUE) {
  assay_imputed <- SummarizedExperiment::assay(se_obj)
  for (i in seq_along(batch_groups)) {
    batch_name    <- names(batch_groups)[i]
    batch_samples <- batch_groups[[i]]
    batch_imp <- mp_impute_single_assay(
      se_obj[, batch_samples, drop = FALSE], method,
      paste0("Batch '", batch_name, "'"), seed = seeds[i],
      impute_with_value = impute_with_value, verbose = verbose)
    assay_imputed[, batch_samples] <- SummarizedExperiment::assay(batch_imp)
    mp_log(verbose, "  Batch '%s': imputation complete for %d sample(s).",
           batch_name, length(batch_samples))
  }
  SummarizedExperiment::assay(se_obj) <- assay_imputed
  se_obj
}

# --- transformation, normalisation, batch correction -----------------------

mp_transform_and_normalise <- function(q, assay_name, normalisation_method,
                                       plot_qc_charts, legend_mode, verbose) {
  plots <- list()
  if (length(normalisation_method) == 1) {
    mp_log(verbose, "Applying normalisation '%s'.", normalisation_method[1])
    q <- mp_run_quiet(
      QFeatures::normalize(q, i = assay_name, name = "normalised",
                           method = normalisation_method[1]),
      suppress_warnings = FALSE)
  } else {
    mp_log(verbose, "Applying transform '%s' then normalisation '%s'.",
           normalisation_method[1], normalisation_method[2])
    x <- SummarizedExperiment::assay(q[[assay_name]])
    mp_assert_log_domain(x, normalisation_method[1])
    x <- switch(normalisation_method[1],
                log2 = log2(x), log = log(x), log10 = log10(x))
    mp_assert_finite(x, paste("after the", normalisation_method[1], "transform"))
    if (plot_qc_charts) {
      plots$transformed_quantities <-
        mp_plot_density(x, "After transformation", legend_mode)
    }
    transformed_se <- q[[assay_name]]
    SummarizedExperiment::assay(transformed_se) <- x
    q_norm <- QFeatures::QFeatures(list(transformed = transformed_se))
    q_norm <- mp_run_quiet(
      QFeatures::normalize(q_norm, i = "transformed", name = "normalised",
                           method = normalisation_method[2]),
      suppress_warnings = FALSE)
    q[["normalised"]] <- q_norm[["normalised"]]
  }
  mp_assert_finite(SummarizedExperiment::assay(q[["normalised"]]),
                   "after normalisation")
  if (plot_qc_charts) {
    plots$normalised_quantities <- mp_plot_density(
      SummarizedExperiment::assay(q[["normalised"]]),
      "After normalisation", legend_mode)
  }
  list(q = q, assay_name = "normalised", plots = plots)
}

mp_run_combat <- function(x, sample_info, batch_col, covariate_cols, verbose) {
  bio <- mp_align_metadata(sample_info, colnames(x), c(covariate_cols, batch_col))
  mp_assert_finite(x, "before batch correction")
  if (anyNA(x)) {
    stop("Batch correction requires a complete matrix, but ", sum(is.na(x)),
         " value(s) are still missing. Choose an `impute_method` other than ",
         "'none', or set `batch_correction = FALSE`.", call. = FALSE)
  }

  batch <- factor(as.character(bio[[batch_col]]))
  singletons <- names(which(table(batch) < 2))
  if (length(singletons) > 0) {
    # ComBat quietly switches to mean-only correction here and says so with a
    # message that used to be suppressed along with its other output.
    mp_log(verbose,
           "  Note: batch(es) %s contain a single sample. ComBat will correct ",
           paste(singletons, collapse = ", "))
    mp_log(verbose,
           "  means only, not variances, for the whole dataset.")
  }

  mod <- NULL
  if (!is.null(covariate_cols)) {
    if (anyNA(bio[, covariate_cols, drop = FALSE])) {
      stop("`batch_correction_variable_columns` contains missing values for ",
           "one or more assay samples.", call. = FALSE)
    }
    mod <- stats::model.matrix(stats::reformulate(covariate_cols), data = bio)
  }

  corrected <- mp_run_quiet(
    sva::ComBat(dat = x, batch = batch, mod = mod, par.prior = TRUE,
                prior.plots = FALSE),
    suppress_warnings = FALSE)

  stopifnot(identical(dim(corrected), dim(x)),
            identical(colnames(corrected), colnames(x)),
            identical(rownames(corrected), rownames(x)))
  mp_assert_finite(corrected, "after batch correction")
  corrected
}

# --- QC plots --------------------------------------------------------------

mp_plot_density <- function(x, title, legend_mode) {
  show_legend <- switch(legend_mode,
                        show = TRUE,
                        hide = FALSE,
                        auto = ncol(x) <= 12)
  limma::plotDensities(x, main = title,
                       legend = if (show_legend) "topright" else FALSE)
  grDevices::recordPlot()
}

# Joins PCA scores to metadata on the explicit Sample column. The previous
# version rebuilt the key from row names, which tibble metadata does not
# carry, so both PCA plots silently disappeared for anyone using read_csv().
mp_build_batch_pca_plot <- function(assay_matrix, sample_info, batch_col,
                                    title, context_label, verbose) {
  x_complete <- assay_matrix[stats::complete.cases(assay_matrix), , drop = FALSE]
  if (nrow(x_complete) > 0) {
    keep <- apply(x_complete, 1,
                  function(v) all(is.finite(v)) && stats::var(v) > 0)
    x_complete <- x_complete[keep, , drop = FALSE]
  }
  if (nrow(x_complete) < 2) {
    mp_log(verbose, "Skipping %s PCA: fewer than 2 non-constant complete features remain.",
           context_label)
    return(NULL)
  }

  x_pca <- t(x_complete)
  max_rank <- min(10L, ncol(x_pca), nrow(x_pca) - 1L)
  if (max_rank < 2L) {
    mp_log(verbose, "Skipping %s PCA: need at least 2 components (samples=%d, features=%d).",
           context_label, nrow(x_pca), ncol(x_pca))
    return(NULL)
  }

  fit <- stats::prcomp(x_pca, scale. = TRUE, center = TRUE, rank. = max_rank)
  var_pct <- 100 * fit$sdev^2 / sum(fit$sdev^2)

  scores <- data.frame(Sample = rownames(fit$x),
                       PC1 = fit$x[, 1], PC2 = fit$x[, 2],
                       row.names = NULL, stringsAsFactors = FALSE)
  bio <- mp_align_metadata(sample_info, scores$Sample, batch_col)
  scores$Batch <- factor(as.character(bio[[batch_col]]))

  ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, colour = Batch)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(title = title, colour = batch_col,
                  x = sprintf("PC1 (%.1f%%)", var_pct[1]),
                  y = sprintf("PC2 (%.1f%%)", var_pct[2])) +
    ggplot2::theme_classic()
}
