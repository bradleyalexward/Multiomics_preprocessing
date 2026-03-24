# Multiomics.Preprocessing

`Multiomics.Preprocessing` is an R package for preprocessing proteomics and metabolomics data
for downstream multi-omics analysis.

It currently provides two main workflows:

- `prepare_proteins()` for precursor-to-protein preprocessing
- `prepare_metabolites()` for peak-to-metabolite preprocessing

The package supports:

- missing-value handling and imputation
- MNAR or MAR classification for mixed imputation
- normalisation and optional transformation
- optional batch-wise imputation and ComBat batch correction
- feature aggregation and QC plotting

## Installation

```r
remotes::install_github("bradleyalexward/Multiomics_preprocessing")
```

## Development checks

From the package root:

```r
roxygen2::roxygenise()
testthat::test_dir("tests/testthat")
```

## Repository layout

- `R/` contains the package functions
- `tests/testthat/` contains package tests
- `inst/extdata/` contains example input files
- `inst/scripts/` contains helper scripts for broader manual testing
