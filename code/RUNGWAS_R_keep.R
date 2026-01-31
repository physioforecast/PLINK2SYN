#!/usr/bin/env Rscript

# =====================================================
# R GWAS pipeline mimicking PLINK 2.0 --glm
# Memory-efficient, batch-enabled, Firth regression safe
# Supports optional keep file (#IID)
# =====================================================

# ------------------- Dependencies -------------------
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

install_if_missing("data.table")
install_if_missing("logistf")

library(data.table)
library(logistf)

# ------------------- Prepare data -------------------
prepare_data <- function(effects_file,
                         pheno_file,
                         effects_start_col = 1,
                         iid_col = "IID",
                         pheno_col,
                         covar_cols = NULL,
                         keep_file = NULL,        # <-- new argument
                         remove_hash_FID = TRUE,
                         pheno_tab_delim = TRUE) {
  
  # --- Load phenotype/covariates ---
  sep <- ifelse(pheno_tab_delim, "\t", ",")
  Y_all <- fread(pheno_file, sep = sep, data.table = FALSE)
  colnames(Y_all) <- trimws(colnames(Y_all))
  if (remove_hash_FID) colnames(Y_all) <- gsub("^#", "", colnames(Y_all))
  Y_all[[iid_col]] <- as.character(Y_all[[iid_col]])
  
  stopifnot(iid_col %in% colnames(Y_all))
  stopifnot(pheno_col %in% colnames(Y_all))
  if (!is.null(covar_cols)) stopifnot(all(covar_cols %in% colnames(Y_all)))
  
  merged_Y <- Y_all[, c(iid_col, pheno_col, covar_cols), drop=FALSE]
  
  # --- Apply keep file filter if provided ---
  if (!is.null(keep_file)) {
    cat(sprintf("📂 Applying keep file: %s\n", keep_file))
    keep_df <- fread(keep_file, header = TRUE, data.table = FALSE)
    
    # Ensure #IID column present
    if (!("#IID" %in% colnames(keep_df))) {
      stop("❌ Keep file must have header '#IID'")
    }
    keep_iids <- as.character(keep_df[["#IID"]])
    
    before_n <- nrow(merged_Y)
    merged_Y <- merged_Y[merged_Y[[iid_col]] %in% keep_iids, , drop = FALSE]
    after_n <- nrow(merged_Y)
    
    cat(sprintf("✅ Filtered by keep file: %d → %d individuals retained\n",
                before_n, after_n))
  }
  
  # --- Load SNP headers only ---
  X_all <- fread(effects_file, nrows = 0)
  colnames(X_all) <- trimws(colnames(X_all))
  predictor_cols <- colnames(X_all)[effects_start_col:ncol(X_all)]
  
  cat(sprintf("✅ Loaded phenotype/covariates: %d samples\n", nrow(merged_Y)))
  cat(sprintf("✅ Detected %d SNPs in effects file\n", length(predictor_cols)))
  
  return(list(merged_Y = merged_Y, snp_cols = predictor_cols))
}

# ------------------- GWAS regression -------------------
run_gwas <- function(effects_file,
                     pheno_file,
                     effects_start_col = 7,
                     iid_col = "IID",
                     pheno_col,
                     covar_cols = NULL,
                     keep_file = NULL,                # <-- new argument
                     regression_type = c("linear", "logistic", "firth"),
                     batch_size = 500,
                     output_prefix = "results/test_run") {
  
  regression_type <- match.arg(regression_type)
  
  # Prepare data
  data_prep <- prepare_data(effects_file, pheno_file,
                            effects_start_col, iid_col, pheno_col, covar_cols,
                            keep_file = keep_file)
  merged_Y <- data_prep$merged_Y
  snp_cols <- data_prep$snp_cols
  n_snps <- length(snp_cols)
  
  results <- data.frame(SNP = character(),
                        beta = numeric(),
                        se = numeric(),
                        p = numeric(),
                        stringsAsFactors = FALSE)
  
  # Split SNPs into batches
  if (is.null(batch_size)) {
    batches <- list(snp_cols)
  } else {
    batches <- split(snp_cols, ceiling(seq_along(snp_cols)/batch_size))
  }
  
  snp_counter <- 0
  for (batch in batches) {
    snp_df <- tryCatch({
      fread(effects_file, select = c(iid_col, batch), data.table = FALSE)
    }, error = function(e) {
      warning(sprintf("Skipping batch %s: %s", paste(batch, collapse=", "), e$message))
      return(NULL)
    })
    if (is.null(snp_df)) next
    snp_df[[iid_col]] <- as.character(snp_df[[iid_col]])
    
    # Merge with phenotype/covariates
    merged_batch <- merge(merged_Y, snp_df, by = iid_col)
    
    for (snp in batch) {
      snp_counter <- snp_counter + 1
      # cat(sprintf("Processing SNP %d / %d (%.2f%%): %s\n", 
      #             snp_counter, n_snps, 100*snp_counter/n_snps, snp))
      # flush.console()
      
      y_sub <- merged_batch[[pheno_col]]
      x_snp <- merged_batch[[snp]]
      X_covars_sub <- if (!is.null(covar_cols)) as.data.frame(merged_batch[, covar_cols, drop=FALSE]) else NULL
      
      # Remove missing
      covar_na <- if (!is.null(X_covars_sub)) apply(is.na(X_covars_sub), 1, any) else rep(FALSE, length(y_sub))
      keep <- !is.na(y_sub) & !is.na(x_snp) & !covar_na
      y_sub <- y_sub[keep]
      x_snp <- x_snp[keep]
      if (!is.null(X_covars_sub)) X_covars_sub <- X_covars_sub[keep, , drop=FALSE]
      if (length(y_sub) < 10) next
      
      # Recode phenotype for logistic/firth
      if (regression_type %in% c("logistic", "firth")) {
        unique_vals <- sort(unique(y_sub))
        if (!all(unique_vals %in% c(0, 1))) {
          max_val <- max(unique_vals)
          y_sub <- ifelse(y_sub == max_val, 1, 0)
        }
      }
      
      df <- data.frame(y = y_sub, SNP = x_snp)
      if (!is.null(X_covars_sub)) df <- cbind(df, X_covars_sub)
      
      beta <- se <- pval <- NA
      tryCatch({
        if (regression_type == "linear") {
          fit <- lm(y ~ ., data = df)
          s <- summary(fit)
          if ("SNP" %in% rownames(s$coefficients)) {
            beta <- s$coefficients["SNP", 1]
            se   <- s$coefficients["SNP", 2]
            pval <- s$coefficients["SNP", 4]
          }
        } else if (regression_type == "logistic") {
          fit <- glm(y ~ ., family = binomial(), data = df)
          s <- summary(fit)
          if ("SNP" %in% rownames(s$coefficients)) {
            beta <- s$coefficients["SNP", 1]
            se   <- s$coefficients["SNP", 2]
            pval <- s$coefficients["SNP", 4]
          }
        } else if (regression_type == "firth") {
          formula_str <- paste("y ~ SNP", 
                               if (!is.null(covar_cols)) paste("+", paste(covar_cols, collapse=" + ")) else "")
          fit <- logistf::logistf(as.formula(formula_str), data = df, pl=TRUE)
          beta <- fit$coefficients["SNP"]
          se <- sqrt(diag(fit$var))["SNP"]
          if (is.na(se)) se <- (fit$ci.upper["SNP"] - fit$ci.lower["SNP"]) / (2 * 1.96)
          pval <- fit$prob["SNP"]
        }
      }, error = function(e) {
        cat("\n❌ Regression failed for SNP:", snp, "\n")
        cat("   Error message:", e$message, "\n")
      })
      
      results <- rbind(results, data.frame(SNP = snp, beta = beta, se = se, p = pval))
    }
  }
  
  # --- Save results ---
  out_dir <- dirname(output_prefix)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_csv <- paste0(output_prefix, "_", regression_type, ".csv")
  write.csv(results, out_csv, row.names = FALSE)
  cat("✅ Saved CSV: ", out_csv, "\n")
  cat("✅ All SNPs processed.\n")
  
  return(results)
}

# ------------------- Example usage -------------------
# Adjust paths and filenames to your system
# path <- Sys.getenv("testdatapath")     # From config.sh
# subpath <- Sys.getenv("resultsdir")    # From config.sh
# fname <- "1kgp3_50k_yesmiss_Av_nonintdose"

# n <- 1000
# effects_file <- paste0(path, fname, "_recode_varIDs_A.raw")
# pheno_file   <- paste0(path, subpath,fname, "_combined_phenocov.csv")
# keep_file    <- paste0(path, subpath,fname, "_recode_varIDs_subset_",n,".keep")  # must have #IID header

# model <- "firth"
# phenotype <- "ybool"
# cov_ <- c("COV_1")

# results <- run_gwas(
#   effects_file      = effects_file,
#   pheno_file        = pheno_file,
#   keep_file         = keep_file,          # <-- subset of individuals
#   effects_start_col = 7,
#   iid_col           = "IID",
#   pheno_col         = phenotype,
#   covar_cols        = cov_,
#   regression_type   = model,
#   batch_size        = 20000,
#   output_prefix     = paste0(path, sbpath,fname, "_", phenotype, "_", paste(cov_, collapse="_"), "_glm_", model, "_keep",n)
# )


# ------------------- Setup -------------------
# Read paths from environment variables (set in config.sh)
path <- Sys.getenv("testdatapath")
sbpath <- Sys.getenv("resultsdir")

# Validate required variables
if (path == "" || sbpath == "") {
  stop("ERROR: testdatapath and resultsdir must be set as environment variables.\n",
       "Make sure config.sh is sourced before running this script.\n",
       "Example in config.sh:\n",
       "  export testdatapath='test_data/'\n",
       "  export resultsdir='PLINK2SYN_GLM_TESTS_YYYYMMDD/'")
}

fnames <- c(
  "1kgp3_50k_nomiss_Av_nonintdose",
  "1kgp3_50k_yesmiss_Av_nonintdose"
)

# Phenotypes
phenotypes <- c("y", "ybool")

# Subset sizes / permutations
thin <- c(1000, 32017, 32000)

# Covariate options
cov_options <- list(
  c("COV_1"),  # with covariate
  NULL         # no covariates
)

batch_size <- 20000


fnames <- c(
  "1kgp3_50k_nomiss_Av_nonintdose",
  "1kgp3_50k_yesmiss_Av_nonintdose"
)

phenotypes <- c("y", "ybool")
thin <- c(1000, 32017, 32000)
cov_options <- list(
  c("COV_1"),  # with covariate
  NULL         # no covariates
)
batch_size <- 20000

# ------------------- Full loop with early exit after 1 iteration -------------------
done <- FALSE  # flag to stop after first iteration

for (fname in fnames) {
  if (done) break
  for (n in thin) {
    if (done) break
    for (cov_ in cov_options) {
      if (done) break
      for (phenotype in phenotypes) {
        if (done) break
        
        # Determine regression types based on phenotype
        regressions <- if (phenotype == "y") "linear" else c("logistic", "firth")
        
        for (model in regressions) {
          
          effects_file <- paste0(path, fname, "_recode_varIDs_A.raw")
          pheno_file   <- paste0(path, sbpath, fname, "_combined_phenocov.csv")
          keep_file    <- paste0(path, sbpath, fname, "_recode_varIDs_subset_", n, ".keep")
          
          cov_suffix <- if (is.null(cov_)) "noCov" else paste(cov_, collapse="_")
          
          output_prefix <- paste0(path, sbpath, fname, "_", phenotype, "_", cov_suffix,
                                  "_glm_", model, "_keep", n)
          
          cat(sprintf("\n📌 Running GWAS for fname=%s, n=%d, cov=%s, phenotype=%s, model=%s\n", 
                      fname, n, cov_suffix, phenotype, model))
          
          results <- run_gwas(
            effects_file      = effects_file,
            pheno_file        = pheno_file,
            keep_file         = keep_file,
            effects_start_col = 7,
            iid_col           = "IID",
            pheno_col         = phenotype,
            covar_cols        = cov_,
            regression_type   = model,
            batch_size        = batch_size,
            output_prefix     = output_prefix
          )
          
          cat(sprintf("✅ Completed GWAS for fname=%s, n=%d, cov=%s, phenotype=%s, model=%s\n", 
                      fname, n, cov_suffix, phenotype, model))
          
          # done <- TRUE  # stop after this first iteration
          # break       # break the innermost loop
        }
      }
    }
  }
}
