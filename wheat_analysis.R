# Wheat genomic prediction across environments
# Run this script from the RStudio project root.

# Required package
if (!requireNamespace("rrBLUP", quietly = TRUE)) {
  stop('Install the required package first with install.packages("rrBLUP").')
}

# Input and output folders
data_dir <- "data"
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

required_files <- c(
  "GenoInfo.csv",
  "PhenoInfo_TPE1.csv",
  "PhenoInfo_TPE2.csv",
  "PhenoInfo_TPE3.csv"
)
missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
if (length(missing_files) > 0) {
  stop(
    "Missing input file(s) in data/: ",
    paste(missing_files, collapse = ", "),
    ". See data/README.md."
  )
}

# Load phenotype and marker data
tpe1 <- read.csv(file.path(data_dir, "PhenoInfo_TPE1.csv"), stringsAsFactors = FALSE)
tpe2 <- read.csv(file.path(data_dir, "PhenoInfo_TPE2.csv"), stringsAsFactors = FALSE)
tpe3 <- read.csv(file.path(data_dir, "PhenoInfo_TPE3.csv"), stringsAsFactors = FALSE)
geno <- read.csv(
  file.path(data_dir, "GenoInfo.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

tpe1$TPE <- "TPE1"
tpe2$TPE <- "TPE2"
tpe3$TPE <- "TPE3"
pheno_all <- rbind(tpe1, tpe2, tpe3)

# Basic data audit
cat("Rows per TPE:\n")
print(c(TPE1 = nrow(tpe1), TPE2 = nrow(tpe2), TPE3 = nrow(tpe3)))
cat("Missing values by TPE:\n")
print(list(TPE1 = colSums(is.na(tpe1)), TPE2 = colSums(is.na(tpe2)), TPE3 = colSums(is.na(tpe3))))
cat("Records per site:\n")
print(list(TPE1 = table(tpe1$Site), TPE2 = table(tpe2$Site), TPE3 = table(tpe3$Site)))
cat("Duplicate Site/line records:\n")
print(c(
  TPE1 = sum(duplicated(tpe1[c("Site", "CIMMYT_ID")])),
  TPE2 = sum(duplicated(tpe2[c("Site", "CIMMYT_ID")])),
  TPE3 = sum(duplicated(tpe3[c("Site", "CIMMYT_ID")]))
))
phenotype_ids <- unique(pheno_all$CIMMYT_ID)
cat("Phenotype line IDs found in marker data:", sum(phenotype_ids %in% geno$GID), "of", length(phenotype_ids), "\n")

# Build the genomic relationship matrix (G)
marker_matrix <- 1 - as.matrix(geno[, -1, drop = FALSE])
storage.mode(marker_matrix) <- "numeric"
allele_frequency <- colMeans(marker_matrix) / 2
keep_markers <- allele_frequency > 0 & allele_frequency < 1
marker_matrix <- marker_matrix[, keep_markers, drop = FALSE]
allele_frequency <- allele_frequency[keep_markers]
marker_centered <- scale(marker_matrix, center = TRUE, scale = FALSE)
denominator <- sum(2 * allele_frequency * (1 - allele_frequency))
G <- tcrossprod(marker_centered) / denominator
rownames(G) <- geno$GID
colnames(G) <- geno$GID

# Cross-validation helper
run_gblup_cv <- function(dat, G, scenario = c("CV1", "CV2"), seed = 1, folds = 5) {
  scenario <- match.arg(scenario)
  set.seed(seed)

  site_levels <- sort(unique(dat$Site))
  id_levels <- rownames(G)
  ids <- unique(dat$CIMMYT_ID)

  if (scenario == "CV1") {
    id_fold <- sample(rep(seq_len(folds), length.out = length(ids)))
    names(id_fold) <- ids
    fold_id <- unname(id_fold[dat$CIMMYT_ID])
  } else {
    fold_id <- integer(nrow(dat))
    for (id in ids) {
      rows <- which(dat$CIMMYT_ID == id)
      fold_id[rows] <- sample(rep(seq_len(folds), length.out = length(rows)))
    }
  }

  results <- data.frame(
    fold = seq_len(folds),
    site_mean_RMSE = NA_real_,
    GBLUP_RMSE = NA_real_
  )

  for (k in seq_len(folds)) {
    train <- dat[fold_id != k, ]
    test <- dat[fold_id == k, ]

    site_means <- aggregate(YLD ~ Site, data = train, mean)
    pred_site <- site_means$YLD[match(test$Site, site_means$Site)]

    train$Site <- factor(train$Site, levels = site_levels)
    test$Site <- factor(test$Site, levels = site_levels)
    X_train <- model.matrix(~ Site, data = train)
    X_test <- model.matrix(~ Site, data = test)
    Z_train <- model.matrix(~ 0 + factor(CIMMYT_ID, levels = id_levels), data = train)
    Z_test <- model.matrix(~ 0 + factor(CIMMYT_ID, levels = id_levels), data = test)

    fit <- rrBLUP::mixed.solve(y = train$YLD, X = X_train, Z = Z_train, K = G)
    pred_gblup <- as.vector(X_test %*% fit$beta + Z_test %*% fit$u)

    results$site_mean_RMSE[k] <- sqrt(mean((test$YLD - pred_site)^2))
    results$GBLUP_RMSE[k] <- sqrt(mean((test$YLD - pred_gblup)^2))
  }

  results
}

# Repeat each validation scenario with two fixed random seeds
cv_runs <- list(
  TPE1_CV1 = list(run_gblup_cv(tpe1, G, "CV1", seed = 789), run_gblup_cv(tpe1, G, "CV1", seed = 101)),
  TPE1_CV2 = list(run_gblup_cv(tpe1, G, "CV2", seed = 123), run_gblup_cv(tpe1, G, "CV2", seed = 456)),
  TPE2_CV1 = list(run_gblup_cv(tpe2, G, "CV1", seed = 202), run_gblup_cv(tpe2, G, "CV1", seed = 408)),
  TPE2_CV2 = list(run_gblup_cv(tpe2, G, "CV2", seed = 204), run_gblup_cv(tpe2, G, "CV2", seed = 405)),
  TPE3_CV1 = list(run_gblup_cv(tpe3, G, "CV1", seed = 301), run_gblup_cv(tpe3, G, "CV1", seed = 406)),
  TPE3_CV2 = list(run_gblup_cv(tpe3, G, "CV2", seed = 302, folds = 3), run_gblup_cv(tpe3, G, "CV2", seed = 407, folds = 3))
)

summary_results <- do.call(rbind, lapply(names(cv_runs), function(nm) {
  both_runs <- do.call(rbind, cv_runs[[nm]])
  baseline <- mean(both_runs$site_mean_RMSE)
  gblup <- mean(both_runs$GBLUP_RMSE)
  data.frame(
    TPE = substr(nm, 1, 4),
    Scenario = if (grepl("CV1", nm)) "CV1: entirely new lines" else "CV2: known lines, missing site record",
    Baseline_RMSE = baseline,
    GBLUP_RMSE = gblup,
    RMSE_reduction_pct = round(100 * (baseline - gblup) / baseline, 1)
  )
}))
rownames(summary_results) <- NULL
print(summary_results)
write.csv(summary_results, "results/summary_results.csv", row.names = FALSE)

# Plot the updated summary
cv2_plot <- summary_results[grepl("^CV2", summary_results$Scenario), ]
cv1_plot <- summary_results[grepl("^CV1", summary_results$Scenario), ]
cv2_plot <- cv2_plot[order(cv2_plot$TPE), ]
cv1_plot <- cv1_plot[order(cv1_plot$TPE), ]

png("figures/prediction_comparison.png", width = 1800, height = 900, res = 150)
old_par <- par(mfrow = c(1, 2), mar = c(5, 4, 4, 1))
barplot(
  cv2_plot$RMSE_reduction_pct,
  names.arg = cv2_plot$TPE,
  ylim = c(-5, 35),
  main = "Known lines: CV2",
  ylab = "% RMSE reduction",
  col = "steelblue"
)
abline(h = 0)
barplot(
  cv1_plot$RMSE_reduction_pct,
  names.arg = cv1_plot$TPE,
  ylim = c(-6, 5),
  main = "New lines: CV1",
  ylab = "% RMSE reduction",
  col = "darkseagreen3"
)
abline(h = 0)
par(old_par)
dev.off()
