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
#' @param zero_handling Either `"na"` (replace 0 with NA) or `"leave"`.
#' Default is "leave"
#' @param impute_method Missing-value imputation strategy.
#' Can be any of: "none", "mixed", "bpca", "knn", "QRILC", "MLE", "MLE2", "MinDet", "MinProb", "min", "zero", "nbavg", "with", "RF".
#' When mixed is selected, MNAR features imputed via QRILC, MAR features imputed via RF.
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
#' @param impute_behaviour Whether to impute missing values globally across all samples (`"global"`)
#' or separately within each batch using `batch_column` (`"per batch"`).
#' Default is `"global"`. For `impute_method = "mixed"`, retained features must have at least one observed value in each batch when `"per batch"` is used.
#' @param missing_feature_proportion_threshold Length 1 or 2 numeric threshold(s) in `[0, 1]`.
#' Default is c(0.3, 0.75). Here, features with more than 30% missingness are removed (first numeric) unless they are MNAR, in which case they are retained up to 75% missingness.
#' If only one numeric is supplied like c(0.3), this will be the missingness threshold of all features.
#' Two numbers must only be supplied if impute_method = "mixed"
#' @param normalisation_method Length 1 normalisation method, or length 2 transform+method. Transform methods include log, log2, log10.
#' Normalise methods include "sum", "max", "center.mean", "center.median", "div.mean", "div.median", "diff.meda", "quantiles", "quantiles.robust" or "vsn".
#' Example: c("vsn) or c("Log10, "center.median").
#' Default is vsn.
#' @param batch_correction Apply ComBat batch correction.
#' @param batch_column Column in `sample_meta_data` identifying which sample belongs to which batch.
#' Used for optional batch-wise imputation and ComBat batch correction.
#' @param batch_correction_variable_columns Optional biological covariates of interest to preserve during batch correction
#' @param protein_aggregator_method Protein aggregation function. Function by which peptides are aggregated to protein level.
#' Default is MsCoreUtils::robustSummary, other examples include MsCoreUtils::medianPolish(), base::colMeans(), matrixStats::colMedians(), etc.
#' Function must take a matrix as input and return a vector of length equal to the column length of the input.
#' @param protein_aggregator_column Column in precursor-level data.frames mapping each precursor to thier identified protein (must share same name).
#' @param plot_qc_charts Return QC plots included density plots of raw, normalised, and transformed abundances, and PCA before and after batch correction.
#' @param show_plot_legends Legend behavior for density plots (`"auto"`, `"show"`, `"hide"`).
#' Default is auto, which hides legends if more than 12 samples are present.
#' @param verbose If `TRUE`, print step-by-step progress messages.
#' Default is TRUE.
#'
#' @return Protein abundance data.frame, or a list with proteins plus optional extras.
#' @export
prepare_proteins <- function(input_files,
                             fraction_file_names = NULL,
                             precursor_abundance_columns,
                             precursor_id_column = "Precursor.Id",
                             sample_meta_data,
                             zero_handling = c("leave", "na"),
                             impute_method = c("none", "mixed", "bpca", "knn", "QRILC", "MLE", "MLE2", "MinDet", "MinProb", "min", "zero", "nbavg", "with", "RF"),
                             mnar_variables = NULL,
                             mnar_significance_threshold = 0.01,
                             return_mnar_results = TRUE,
                             impute_with_value = NULL,
                             impute_behaviour = c("global", "per batch"),
                             missing_feature_proportion_threshold = c(0.3, 0.75),
                             normalisation_method = c("vsn"),
                             batch_correction = FALSE,
                             batch_column = "Batch",
                             batch_correction_variable_columns = NULL,
                             protein_aggregator_method = MsCoreUtils::robustSummary,
                             protein_aggregator_column = "Protein.Group",
                             plot_qc_charts = FALSE,
                             show_plot_legends = c("auto", "show", "hide"),
                             verbose = TRUE) {

  zero_handling <- match.arg(zero_handling)
  show_plot_legends <- match.arg(show_plot_legends)
  impute_method <- match.arg(impute_method)
  impute_behaviour <- match.arg(impute_behaviour)

  log_step <- function(fmt, ...) {
    if (isTRUE(verbose)) {
      message(sprintf(fmt, ...))
    }
  }

  run_quiet <- function(expr, suppress_warnings = TRUE) {
    value <- NULL
    if (suppress_warnings) {
      suppressWarnings(
        suppressMessages(
          utils::capture.output(value <- eval.parent(substitute(expr)), type = "output")
        )
      )
    } else {
      suppressMessages(
        utils::capture.output(value <- eval.parent(substitute(expr)), type = "output")
      )
    }
    value
  }

  resolve_quant_cols <- function(df, cols, fraction_label) {
    if (is.numeric(cols)) {
      idx <- unique(as.integer(cols))
      if (any(is.na(idx)) || any(idx < 1 | idx > ncol(df))) {
        stop("`precursor_abundance_columns` contains invalid numeric indices for fraction '", fraction_label, "'.")
      }
      return(idx)
    }

    if (is.character(cols)) {
      missing_cols <- setdiff(cols, colnames(df))
      if (length(missing_cols) > 0) {
        stop("`precursor_abundance_columns` contains missing column names for fraction '", fraction_label,
             "': ", paste(missing_cols, collapse = ", "))
      }
      return(match(cols, colnames(df)))
    }

    if (is.logical(cols)) {
      if (length(cols) != ncol(df)) {
        stop("Logical `precursor_abundance_columns` must have length ncol(df) for fraction '", fraction_label, "'.")
      }
      idx <- which(cols)
      if (length(idx) == 0) {
        stop("Logical `precursor_abundance_columns` selected no columns for fraction '", fraction_label, "'.")
      }
      return(idx)
    }

    stop("`precursor_abundance_columns` must be numeric, character, or logical.")
  }

  coerce_abundance_to_numeric <- function(df, quant_cols, fraction_label) {
    for (j in quant_cols) {
      if (!is.numeric(df[[j]])) {
        converted <- suppressWarnings(as.numeric(df[[j]]))
        bad <- !is.na(df[[j]]) & is.na(converted)
        if (any(bad)) {
          stop("Non-numeric values found in abundance column '", colnames(df)[j], "' for fraction '",
               fraction_label, "'.")
        }
        df[[j]] <- converted
      }
    }
    df
  }

  detect_sample_col <- function(df, assay_sample_names, label) {
    col_matches <- vapply(df, function(x) sum(as.character(x) %in% assay_sample_names), numeric(1))
    if (all(col_matches == 0)) {
      stop("Could not identify a sample ID column in ", label, " by matching assay sample names.")
    }
    names(which.max(col_matches))
  }

  build_sample_info <- function(df, assay_sample_names, label) {
    sample_col <- detect_sample_col(df, assay_sample_names, label)
    meta <- df
    meta$Sample <- as.character(meta[[sample_col]])

    if (anyNA(meta$Sample) || any(meta$Sample == "")) {
      stop(label, " contains missing or empty sample identifiers.")
    }
    if (anyDuplicated(meta$Sample) > 0) {
      stop(label, " contains duplicated sample identifiers in column '", sample_col, "'.")
    }

    missing_samples <- setdiff(assay_sample_names, meta$Sample)
    if (length(missing_samples) > 0) {
      stop(label, " is missing sample IDs required by the assay: ",
           paste(missing_samples, collapse = ", "))
    }

    meta <- meta[match(assay_sample_names, meta$Sample), , drop = FALSE]
    rownames(meta) <- meta$Sample
    attr(meta, "sample_col") <- sample_col
    meta
  }

  remove_missing_constant <- function(df, quant_cols) {
    x <- df[, quant_cols, drop = FALSE]
    all_missing <- rowSums(is.na(x)) == ncol(x)
    constant <- apply(x, 1, function(v) {
      v <- v[!is.na(v)]
      length(v) >= 2 && stats::var(v) == 0
    })
    keep <- !(all_missing | constant)
    list(
      data = df[keep, , drop = FALSE],
      n_missing = sum(all_missing),
      n_constant = sum(constant),
      n_before = nrow(df),
      n_after = sum(keep)
    )
  }

  plot_density <- function(x, title, legend_mode) {
    show_legend <- switch(legend_mode,
                          show = TRUE,
                          hide = FALSE,
                          auto = ncol(x) <= 12)
    limma::plotDensities(x,
                         main = title,
                         legend = if (show_legend) "topright" else FALSE)
    grDevices::recordPlot()
  }

  build_batch_pca_plot <- function(assay_matrix, bio_df, batch_col, title, context_label) {
    x_complete <- assay_matrix[stats::complete.cases(assay_matrix), , drop = FALSE]
    if (nrow(x_complete) < 2) {
      log_step("Skipping %s PCA: fewer than 2 complete features remain after NA filtering.", context_label)
      return(NULL)
    }

    x_pca <- t(x_complete)
    max_rank <- min(10L, ncol(x_pca), nrow(x_pca) - 1L)
    if (max_rank < 2L) {
      log_step("Skipping %s PCA: need at least 2 principal components (samples=%d, complete features=%d).",
               context_label, nrow(x_pca), ncol(x_pca))
      return(NULL)
    }

    pca_fit <- stats::prcomp(x = x_pca, scale. = TRUE, center = TRUE, rank. = max_rank)

    components_pca <- tibble::rownames_to_column(as.data.frame(pca_fit[["x"]]), "Sample")
    components_pca <- dplyr::inner_join(
      components_pca,
      tibble::rownames_to_column(as.data.frame(bio_df), "Sample"),
      by = "Sample"
    )

    if (!all(c("PC1", "PC2") %in% colnames(components_pca))) {
      log_step("Skipping %s PCA: PC1/PC2 were not available from PCA output.", context_label)
      return(NULL)
    }
    if (nrow(components_pca) == 0) {
      log_step("Skipping %s PCA: no overlapping samples between PCA scores and metadata.", context_label)
      return(NULL)
    }

    components_pca$Batch <- as.factor(components_pca[[batch_col]])

    factoextra::fviz_pca_ind(
      pca_fit,
      mean.point = FALSE,
      label = FALSE,
      pointsize = 0
    ) +
      ggplot2::geom_point(
        data = components_pca,
        ggplot2::aes(x = PC1, y = PC2, text = Sample, color = Batch),
        size = 3
      ) +
      ggplot2::labs(title = title, color = batch_col) +
      ggplot2::theme_classic()
  }

  classify_mnar <- function(object,
                            sample_info,
                            covariate_cols,
                            assay_name,
                            feature_name,
                            p_adj_cutoff,
                            mar_threshold) {
    assay_mat <- SummarizedExperiment::assay(object[[assay_name]])
    binary_assay <- !is.na(assay_mat)
    binary_assay <- binary_assay[rowSums(binary_assay) >= 1, , drop = FALSE]

    binary_assay <- as.data.frame(binary_assay)
    binary_assay <- tibble::rownames_to_column(binary_assay, "Feature")
    binary_assay <- tidyr::pivot_longer(binary_assay,
                                        cols = -"Feature",
                                        names_to = "Sample",
                                        values_to = "Identified")

    covariates2 <- as.data.frame(sample_info[, c("Sample", covariate_cols), drop = FALSE])
    binary_assay <- dplyr::left_join(binary_assay, as.data.frame(covariates2), by = "Sample")

    fml <- stats::reformulate(covariate_cols, response = "Identified")

    feats <- unique(binary_assay$Feature)
    pvals <- rep(1, length(feats))
    names(pvals) <- feats
    frac_missing <- rep(NA_real_, length(feats))
    names(frac_missing) <- feats

    for (i in seq_along(feats)) {
      feature_df <- binary_assay[binary_assay$Feature == feats[i], , drop = FALSE]
      fit <- suppressWarnings(try(stats::glm(fml, feature_df, family = "binomial"), silent = TRUE))
      if (!inherits(fit, "try-error")) {
        p <- suppressWarnings(try(stats::coef(summary(fit))[, "Pr(>|z|)"], silent = TRUE))
        if (!inherits(p, "try-error")) {
          p_adj <- suppressWarnings(stats::p.adjust(p, method = "BH"))
          p_adj <- p_adj[is.finite(p_adj)]
          if (length(p_adj) > 0) pvals[i] <- min(p_adj)
        }
      }
      frac_missing[i] <- mean(!feature_df$Identified)
    }

    adj_pvals <- stats::p.adjust(pvals, method = "BH")
    results <- data.frame(
      pvals = pvals,
      adj.pvals = adj_pvals,
      MNAR = adj_pvals <= p_adj_cutoff,
      Frac_missing = frac_missing,
      Remove_feature = (adj_pvals > p_adj_cutoff) & (frac_missing > mar_threshold)
    )
    results <- tibble::rownames_to_column(results, feature_name)

    rowdata_df <- as.data.frame(SummarizedExperiment::rowData(object[[assay_name]]))
    if (!(feature_name %in% colnames(rowdata_df))) {
      rowdata_df[[feature_name]] <- rownames(rowdata_df)
    }
    rowdata_df <- dplyr::left_join(rowdata_df, results, by = feature_name)
    SummarizedExperiment::rowData(object[[assay_name]]) <- rowdata_df

    list(object = object, results = results)
  }

  build_batch_groups <- function(sample_info, assay_sample_names, batch_col) {
    if (!(batch_col %in% colnames(sample_info))) {
      stop("`batch_column` ('", batch_col, "') is not present in `sample_meta_data`.")
    }
    batch_values <- as.character(sample_info[[batch_col]])
    if (anyNA(batch_values) || any(batch_values == "")) {
      stop("`batch_column` contains missing or empty values for one or more assay samples.")
    }
    split(assay_sample_names, batch_values)
  }

  find_batch_all_missing_features <- function(se_obj, batch_groups) {
    x <- SummarizedExperiment::assay(se_obj)
    feature_ids <- rownames(x)
    if (is.null(feature_ids)) {
      feature_ids <- as.character(seq_len(nrow(x)))
    }

    details <- vector("list", length(batch_groups))
    n_details <- 0L

    for (i in seq_along(batch_groups)) {
      batch_name <- names(batch_groups)[i]
      batch_samples <- batch_groups[[i]]
      if (length(batch_samples) == 0) {
        next
      }

      all_missing <- rowSums(!is.na(x[, batch_samples, drop = FALSE])) == 0
      if (any(all_missing)) {
        n_details <- n_details + 1L
        details[[n_details]] <- data.frame(
          Feature = feature_ids[all_missing],
          Batch = batch_name,
          stringsAsFactors = FALSE
        )
      }
    }

    if (n_details == 0L) {
      return(data.frame(Feature = character(0), Batch = character(0), stringsAsFactors = FALSE))
    }

    do.call(rbind, details[seq_len(n_details)])
  }

  impute_single_assay <- function(se_obj, method, label) {
    if (method == "mixed") {
      randna <- !as.logical(SummarizedExperiment::rowData(se_obj)$MNAR)
      randna[is.na(randna)] <- TRUE
      imp <- run_quiet(QFeatures::impute(se_obj, method = "mixed", randna = randna, mar = "RF", mnar = "QRILC"))
    } else if (method == "with") {
      imp <- run_quiet(QFeatures::impute(se_obj, method = "with", val = impute_with_value))
    } else {
      imp <- run_quiet(QFeatures::impute(se_obj, method = method))
    }
    imp
  }

  impute_assay_by_batch <- function(se_obj, batch_groups, method) {
    assay_imputed <- SummarizedExperiment::assay(se_obj)

    for (i in seq_along(batch_groups)) {
      batch_name <- names(batch_groups)[i]
      batch_samples <- batch_groups[[i]]
      batch_se <- se_obj[, batch_samples, drop = FALSE]
      batch_imp <- impute_single_assay(batch_se, method, paste0("Batch '", batch_name, "'"))
      assay_imputed[, batch_samples] <- SummarizedExperiment::assay(batch_imp)
      log_step("  Batch '%s': imputation complete for %d sample(s).", batch_name, length(batch_samples))
    }

    SummarizedExperiment::assay(se_obj) <- assay_imputed
    se_obj
  }

  # -------------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------------
  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }
  if (!is.list(input_files) || inherits(input_files, "data.frame") || length(input_files) == 0) {
    stop("`input_files` must be a non-empty list of one or more data.frames.")
  }
  if (!all(vapply(input_files, inherits, logical(1), what = "data.frame"))) {
    stop("All elements of `input_files` must be data.frames.")
  }
  if (!is.data.frame(sample_meta_data) || nrow(sample_meta_data) == 0) {
    stop("`sample_meta_data` must be a non-empty data.frame.")
  }
  if (!is.character(precursor_id_column) || length(precursor_id_column) != 1 || is.na(precursor_id_column)) {
    stop("`precursor_id_column` must be a single character string.")
  }

  if (!is.null(fraction_file_names)) {
    if (!is.character(fraction_file_names)) stop("`fraction_file_names` must be a character vector when supplied.")
    if (length(fraction_file_names) != length(input_files)) stop("`fraction_file_names` must have one element per input file.")
    if (any(is.na(fraction_file_names) | fraction_file_names == "")) stop("`fraction_file_names` cannot contain NA or empty values.")
    if (any(duplicated(fraction_file_names))) stop("`fraction_file_names` must be unique.")
  }

  if (!is.numeric(mnar_significance_threshold) || length(mnar_significance_threshold) != 1 ||
      is.na(mnar_significance_threshold) || mnar_significance_threshold <= 0 || mnar_significance_threshold > 1) {
    stop("`mnar_significance_threshold` must be a single numeric value in (0, 1].")
  }

  if (!is.numeric(missing_feature_proportion_threshold) ||
      !(length(missing_feature_proportion_threshold) %in% c(1, 2)) ||
      any(is.na(missing_feature_proportion_threshold)) ||
      any(missing_feature_proportion_threshold < 0 | missing_feature_proportion_threshold > 1)) {
    stop("`missing_feature_proportion_threshold` must be numeric length 1 or 2 with values in [0, 1].")
  }
  if (length(missing_feature_proportion_threshold) == 2 &&
      missing_feature_proportion_threshold[2] < missing_feature_proportion_threshold[1]) {
    stop("When length 2, `missing_feature_proportion_threshold` must be c(MAR_threshold, MNAR_threshold) with MNAR >= MAR.")
  }

  if (!is.character(normalisation_method)) stop("`normalisation_method` must be a character vector.")
  normalize_methods <- MsCoreUtils::normalizeMethods()
  if (length(normalisation_method) == 1) {
    if (!(normalisation_method[1] %in% normalize_methods)) {
      stop("Unsupported normalisation method: '", normalisation_method[1], "'.")
    }
  } else if (length(normalisation_method) == 2) {
    if (!(normalisation_method[1] %in% c("log2", "log", "log10"))) {
      stop("When length 2, first normalisation element must be one of: log2, log, log10.")
    }
    if (!(normalisation_method[2] %in% normalize_methods)) {
      stop("Unsupported normalisation method: '", normalisation_method[2], "'.")
    }
  } else {
    stop("`normalisation_method` must have length 1 or 2.")
  }

  if (!is.function(protein_aggregator_method)) stop("`protein_aggregator_method` must be a function.")
  if (!is.character(protein_aggregator_column) || length(protein_aggregator_column) != 1 || is.na(protein_aggregator_column)) {
    stop("`protein_aggregator_column` must be a single character string.")
  }

  if (impute_method == "MLE2") {
    stop("`impute_method = 'MLE2'` is defunct in QFeatures. Please use `impute_method = 'MLE'`.")
  }

  if (impute_method == "mixed") {
    if (is.null(mnar_variables) ||
        !is.character(mnar_variables) ||
        length(mnar_variables) == 0 ||
        any(is.na(mnar_variables) | mnar_variables == "")) {
      stop("`mnar_variables` must be supplied as one or more column names from `sample_meta_data` when `impute_method = 'mixed'`.")
    }
  } else if (!is.null(mnar_variables) &&
             (!is.character(mnar_variables) || any(is.na(mnar_variables) | mnar_variables == ""))) {
    stop("`mnar_variables` must be `NULL` or a character vector of column names from `sample_meta_data`.")
  }

  if (impute_method == "with") {
    if (is.null(impute_with_value) || !is.numeric(impute_with_value) || length(impute_with_value) != 1 || !is.finite(impute_with_value)) {
      stop("When `impute_method = 'with'`, `impute_with_value` must be a single finite numeric value.")
    }
  }

  use_batchwise_imputation <- (impute_method != "none") && identical(impute_behaviour, "per batch")
  needs_batch_info <- use_batchwise_imputation || isTRUE(batch_correction)
  if (needs_batch_info) {
    if (!is.character(batch_column) || length(batch_column) != 1 || is.na(batch_column) || batch_column == "") {
      stop("`batch_column` must be a single non-empty character string when batch-aware steps are enabled.")
    }
    if (!(batch_column %in% colnames(sample_meta_data))) {
      stop("`batch_column` ('", batch_column, "') is not present in `sample_meta_data`.")
    }
  }

  if (batch_correction) {
    if (!is.null(batch_correction_variable_columns)) {
      if (!is.character(batch_correction_variable_columns)) {
        stop("`batch_correction_variable_columns` must be a character vector when supplied.")
      }
      missing_cols <- setdiff(batch_correction_variable_columns, colnames(sample_meta_data))
      if (length(missing_cols) > 0) {
        stop("Missing `batch_correction_variable_columns` in `sample_meta_data`: ", paste(missing_cols, collapse = ", "))
      }
    }
  }

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
    quant_cols_list[[i]] <- resolve_quant_cols(dfs[[i]], precursor_abundance_columns, frac)
    if (precursor_id_column %in% colnames(dfs[[i]])[quant_cols_list[[i]]]) {
      stop("`precursor_id_column` cannot be included in `precursor_abundance_columns` for fraction '", frac, "'.")
    }
    dfs[[i]] <- coerce_abundance_to_numeric(dfs[[i]], quant_cols_list[[i]], frac)
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

  sample_info <- build_sample_info(sample_meta_data, sample_names, "sample_meta_data")
  sample_id_column <- attr(sample_info, "sample_col")
  mnar_variable_columns <- NULL
  if (!is.null(mnar_variables)) {
    mnar_variable_columns <- unique(as.character(mnar_variables))
    missing_mnar_cols <- setdiff(mnar_variable_columns, colnames(sample_info))
    if (length(missing_mnar_cols) > 0) {
      stop("`mnar_variables` contains columns not found in `sample_meta_data`: ",
           paste(missing_mnar_cols, collapse = ", "))
    }
    forbidden_mnar_cols <- intersect(mnar_variable_columns, unique(c("Sample", sample_id_column)))
    if (length(forbidden_mnar_cols) > 0) {
      stop("`mnar_variables` must reference covariate columns in `sample_meta_data`, not the sample ID column: ",
           paste(forbidden_mnar_cols, collapse = ", "))
    }
  }
  batch_groups <- NULL
  if (needs_batch_info) {
    batch_groups <- build_batch_groups(sample_info, sample_names, batch_column)
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

  cleaned <- lapply(seq_along(dfs), function(i) remove_missing_constant(dfs[[i]], quant_cols_list[[i]]))
  for (i in seq_along(cleaned)) {
    dfs[[i]] <- cleaned[[i]]$data
    log_step("  Fraction '%s': removed %d all-missing and %d constant features (kept %d/%d).",
             fraction_file_names[i], cleaned[[i]]$n_missing, cleaned[[i]]$n_constant,
             cleaned[[i]]$n_after, cleaned[[i]]$n_before)
  }

  if (length(dfs) > 1) {
    for (i in seq_along(dfs)) {
      dfs[[i]][[precursor_id_column]] <- paste0(dfs[[i]][[precursor_id_column]], "_", fraction_file_names[i])
    }
  }

  se_list <- list()
  for (i in seq_along(dfs)) {
    se_list[[fraction_file_names[i]]] <- run_quiet(
      QFeatures::readSummarizedExperiment(dfs[[i]],
                                          quantCols = quant_cols_list[[i]],
                                          fnames = precursor_id_column)
    )
  }
  qf <- QFeatures::QFeatures(se_list)
  log_step("Created QFeatures object with assays: %s", paste(names(qf), collapse = ", "))

  # -------------------------------------------------------------------------
  # MNAR/MAR classification and filtering
  # -------------------------------------------------------------------------
  do_mnar <- (impute_method == "mixed")
  mnar_tables <- list()
  initial_assay_names <- names(qf)
  mar_thr <- missing_feature_proportion_threshold[1]
  mnar_thr <- ifelse(length(missing_feature_proportion_threshold) >= 2,
                     missing_feature_proportion_threshold[2], 0.75)

  if (do_mnar) {
    if (is.null(mnar_variable_columns) || length(mnar_variable_columns) == 0) {
      stop("`mnar_variables` must contain at least one covariate column from `sample_meta_data` when `impute_method = 'mixed'`.")
    }
    log_step("Starting MNAR classification for mixed imputation (p_adj <= %.4f).", mnar_significance_threshold)
    for (nm in initial_assay_names) {
      res <- classify_mnar(qf, sample_info, mnar_variable_columns, nm, precursor_id_column, mnar_significance_threshold, mar_thr)
      qf <- res$object
      mnar_tables[[nm]] <- res$results
      log_step("  Assay '%s': MNAR features = %d / %d.", nm, sum(res$results$MNAR, na.rm = TRUE), nrow(res$results))
    }
    if (length(mnar_tables) == 0) {
      stop("MNAR classification returned no results while `impute_method = 'mixed'`.")
    }

    before_counts <- vapply(initial_assay_names, function(nm) nrow(qf[[nm]]), integer(1))
    qf <- run_quiet(QFeatures::filterFeatures(
      qf,
      QFeatures::VariableFilter(field = "Remove_feature", value = 0, condition = "=="),
      i = initial_assay_names
    ))
    qf <- run_quiet(QFeatures::filterFeatures(
      qf,
      QFeatures::VariableFilter(field = "Frac_missing", value = mnar_thr, condition = "<"),
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

  # -------------------------------------------------------------------------
  # Imputation
  # -------------------------------------------------------------------------
  source_assays <- names(qf)
  if (impute_method != "none") {
    log_step("Starting imputation using method '%s'.", impute_method)
    for (nm in source_assays) {
      na_before <- sum(is.na(SummarizedExperiment::assay(qf[[nm]])))
      if (use_batchwise_imputation) {
        if (impute_method == "mixed") {
          batch_all_missing <- find_batch_all_missing_features(qf[[nm]], batch_groups)
          if (nrow(batch_all_missing) > 0) {
            batch_counts <- sort(table(batch_all_missing$Batch), decreasing = TRUE)
            batch_counts_text <- paste(
              sprintf("%s=%d", names(batch_counts), as.integer(batch_counts)),
              collapse = ", "
            )
            example_features <- utils::head(unique(batch_all_missing$Feature), 5)
            stop(
              "Batch-wise mixed imputation cannot be performed because ",
              nrow(batch_all_missing), " feature(s) are completely missing within at least one batch (",
              batch_counts_text, "). Example feature IDs: ",
              paste(example_features, collapse = ", "),
              ". Use a stricter `missing_feature_proportion_threshold` or set `impute_behaviour = 'global'`."
            )
          }
        }
        log_step("  Assay '%s': imputing missing values separately within %d batch(es).", nm, length(batch_groups))
        imp <- impute_assay_by_batch(qf[[nm]], batch_groups, impute_method)
      } else {
        imp <- impute_single_assay(qf[[nm]], impute_method, paste0("Assay '", nm, "'"))
      }

      na_after <- sum(is.na(SummarizedExperiment::assay(imp)))
      if (na_after > 0) {
        stop("Imputation method '", impute_method, "' returned ", na_after,
             " missing values in assay '", nm, "'.")
      }

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
    log_step("Combining %d assays into a single assay.", length(use_assays))
    combined <- do.call(rbind, lapply(use_assays, function(nm) {
      cbind(as.data.frame(SummarizedExperiment::rowData(qf[[nm]])),
            SummarizedExperiment::assay(qf, nm))
    }))
    sample_cols <- colnames(SummarizedExperiment::assay(qf, use_assays[1]))
    se_combined <- run_quiet(QFeatures::readSummarizedExperiment(
      combined,
      quantCols = colnames(combined) %in% sample_cols,
      fnames = precursor_id_column
    ))
    q <- QFeatures::QFeatures(list(se = se_combined))
    se_name <- "se"
  } else {
    q <- QFeatures::QFeatures(list(se = qf[[use_assays[1]]]))
    se_name <- "se"
  }

  # -------------------------------------------------------------------------
  # QC plots
  # -------------------------------------------------------------------------
  plots <- list()
  if (plot_qc_charts) {
    plots$raw_quantities <- plot_density(SummarizedExperiment::assay(q[[se_name]]), "Raw", show_plot_legends)
  }

  # -------------------------------------------------------------------------
  # Transform + normalize
  # -------------------------------------------------------------------------
  if (any(normalisation_method == "vsn")) {
    invisible(vsn::justvsn)
  }
  if (length(normalisation_method) == 1) {
    log_step("Applying normalisation '%s'.", normalisation_method[1])
    q <- run_quiet(QFeatures::normalize(q, i = se_name, name = "norm", method = normalisation_method[1]))
    assay_name <- "norm"
  } else {
    log_step("Applying transform '%s' then normalisation '%s'.", normalisation_method[1], normalisation_method[2])
    x <- SummarizedExperiment::assay(q[[se_name]])
    if (normalisation_method[1] == "log2") x <- log2(x)
    if (normalisation_method[1] == "log") x <- log(x)
    if (normalisation_method[1] == "log10") x <- log10(x)
    if (plot_qc_charts) {
      plots$transformed_quantities <- plot_density(x, "After transformation", show_plot_legends)
    }
    SummarizedExperiment::assay(q[[se_name]]) <- x
    q <- run_quiet(QFeatures::normalize(q, i = se_name, name = "norm", method = normalisation_method[2]))
    assay_name <- "norm"
  }
  if (plot_qc_charts) {
    plots$normalised_quantities <- plot_density(SummarizedExperiment::assay(q[[assay_name]]), "After normalisation", show_plot_legends)
  }

  # -------------------------------------------------------------------------
  # Batch correction
  # -------------------------------------------------------------------------
  if (batch_correction) {
    log_step("Applying ComBat batch correction with batch column '%s'.", batch_column)
    x <- SummarizedExperiment::assay(q[[assay_name]])
    bio <- dplyr::select(sample_info, dplyr::all_of(c(batch_correction_variable_columns, batch_column)))

    if (plot_qc_charts) {
      before_pca_plot <- build_batch_pca_plot(
        assay_matrix = SummarizedExperiment::assay(q[[assay_name]]),
        bio_df = bio,
        batch_col = batch_column,
        title = "Before batch correction",
        context_label = "pre-batch-correction"
      )
      if (!is.null(before_pca_plot)) {
        plots$before_batch_correction <- before_pca_plot
      }
    }

    batch <- as.factor(bio[[batch_column]])
    mod <- NULL
    if (!is.null(batch_correction_variable_columns)) {
      mod <- stats::model.matrix(stats::as.formula(paste("~", paste(batch_correction_variable_columns, collapse = " + "))),
                                 data = bio)
    }
    corrected <- run_quiet(sva::ComBat(dat = x, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE))
    SummarizedExperiment::assay(q[[assay_name]]) <- corrected
    log_step("Batch correction complete.")

    if (plot_qc_charts) {
      after_pca_plot <- build_batch_pca_plot(
        assay_matrix = SummarizedExperiment::assay(q[[assay_name]]),
        bio_df = bio,
        batch_col = batch_column,
        title = "After batch correction",
        context_label = "post-batch-correction"
      )
      if (!is.null(after_pca_plot)) {
        plots$after_batch_correction <- after_pca_plot
      }
    }

  } else {
    log_step("Batch correction skipped.")
  }

  # -------------------------------------------------------------------------
  # Aggregate to proteins
  # -------------------------------------------------------------------------
  row_data <- SummarizedExperiment::rowData(q[[assay_name]])
  if (!(protein_aggregator_column %in% colnames(row_data))) {
    fallback_rowdata <- SummarizedExperiment::rowData(q[[se_name]])
    if (protein_aggregator_column %in% colnames(fallback_rowdata)) {
      row_data[[protein_aggregator_column]] <- fallback_rowdata[[protein_aggregator_column]]
      SummarizedExperiment::rowData(q[[assay_name]]) <- row_data
    } else {
      stop("`protein_aggregator_column` not found in rowData of assay '", assay_name, "'.")
    }
  }

  log_step("Aggregating precursor features to proteins using '%s'.", protein_aggregator_column)
  q <- run_quiet(QFeatures::aggregateFeatures(
    q,
    i = assay_name,
    fcol = protein_aggregator_column,
    name = "Proteins",
    fun = protein_aggregator_method,
    na.rm = TRUE
  ))

  proteins_mat <- SummarizedExperiment::assay(q[["Proteins"]])
  proteins_na <- sum(is.na(proteins_mat))
  if (proteins_na > 0) {
    stop("Protein aggregation produced ", proteins_na, " missing values.")
  }
  proteins <- as.data.frame(proteins_mat)
  log_step("Finished: %d proteins x %d samples.", nrow(proteins), ncol(proteins))

  if (!return_mnar_results && !plot_qc_charts) {
    return(proteins)
  }

  out <- list(proteins = proteins)
  if (return_mnar_results) {
    if (do_mnar && length(mnar_tables) == 0) {
      stop("`mnar_results` is empty while `impute_method = 'mixed'`. MNAR classification failed.")
    }
    out$mnar_results <- mnar_tables
  }
  if (plot_qc_charts) out$plots <- plots
  out
}
