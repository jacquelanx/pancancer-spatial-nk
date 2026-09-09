#!/usr/bin/env Rscript
# ============================================================================
# Cross-platform harmonization and reference-based spatial deconvolution.
#
# Pipeline: per-sample MAD QC -> doublet removal -> within-cancer-type Harmony
#           integration -> RCTD deconvolution -> deconvolution reliability
#           benchmark -> per-cell-type reliability flags.
#
# Three decisions in here determine whether anything downstream is interpretable.
# They are implemented, logged, and written to the run manifest:
#
#   (1) QC thresholds are per-sample and MAD-based. A flat 10% mitochondrial cut
#       deletes real cells in high-metabolic tissue (kidney, liver, heart) and
#       keeps junk in others.
#   (2) Integration is run WITHIN cancer type over patient. In a pan-cancer
#       cohort, batch is confounded with cancer type, so "correcting batch"
#       across types removes the biology the study is about.
#   (3) Deconvolution of a rare cell type is benchmarked, not trusted. NK cells
#       are 0.5-3% of a tumor and a 55um spot holds 1-10 cells; the RCTD weight
#       for NK is not separable from CD8 T at realistic depth. Cell types that
#       fail the simulated-mixture benchmark are flagged unreliable and must not
#       carry a conclusion.
#
# NOTE: this script has not been executed in the environment it was authored in
# (no R toolchain present). Run `Rscript R/harmonize_deconvolve.R --help` and the
# smoke target in the Snakefile before trusting it on a full cohort.
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(harmony)
  library(spacexr)
})

# ------------------------------------------------------------------ logging

log_msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(fmt, ...)))
  flush.console()
}

stop_if <- function(condition, fmt, ...) {
  if (isTRUE(condition)) stop(sprintf(fmt, ...), call. = FALSE)
}

check_versions <- function() {
  v <- as.character(packageVersion("Seurat"))
  log_msg("Seurat %s | SeuratObject %s | spacexr %s", v,
          packageVersion("SeuratObject"), packageVersion("spacexr"))
  stop_if(as.integer(substr(v, 1, 1)) < 5,
          "Seurat v5 required (found %s): the layer API used here does not exist in v4", v)
}

# ------------------------------------------------------------------ QC

#' Per-sample MAD thresholds for a QC metric.
#'
#' Returns the upper (and optionally lower) cut for a metric within one sample,
#' so the threshold adapts to tissue and chemistry instead of being imposed.
mad_thresholds <- function(x, nmads = 3, type = c("higher", "both")) {
  type <- match.arg(type)
  med <- median(x, na.rm = TRUE)
  dev <- mad(x, center = med, na.rm = TRUE)
  if (dev == 0) dev <- .Machine$double.eps  # constant metric: fall back to no filtering
  list(lower = if (type == "both") med - nmads * dev else -Inf,
       upper = med + nmads * dev,
       median = med, mad = dev)
}

qc_filter_sample <- function(obj, sample_id, nmads = 3, min_features = 200) {
  obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent_ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

  mt   <- mad_thresholds(obj$percent_mt, nmads, "higher")
  feat <- mad_thresholds(log10(obj$nFeature_RNA + 1), nmads, "both")
  cnt  <- mad_thresholds(log10(obj$nCount_RNA + 1), nmads, "both")

  keep <- obj$percent_mt <= mt$upper &
    log10(obj$nFeature_RNA + 1) >= feat$lower &
    log10(obj$nFeature_RNA + 1) <= feat$upper &
    log10(obj$nCount_RNA + 1)   >= cnt$lower &
    log10(obj$nCount_RNA + 1)   <= cnt$upper &
    obj$nFeature_RNA >= min_features

  log_msg("QC %s: %d -> %d cells (%.1f%% removed); mt cut %.1f%% (median %.1f%%)",
          sample_id, ncol(obj), sum(keep), 100 * (1 - sum(keep) / ncol(obj)),
          mt$upper, mt$median)

  attr(keep, "thresholds") <- list(
    sample = sample_id, mt_upper = mt$upper, mt_median = mt$median,
    n_before = ncol(obj), n_after = sum(keep)
  )
  keep
}

#' DoubletFinder wrapper tolerant of the package's renamed functions.
#'
#' DoubletFinder renamed paramSweep_v3 -> paramSweep and doubletFinder_v3 ->
#' doubletFinder between releases; pinning is preferable, but an unpinned cluster
#' module should not take the whole run down.
remove_doublets <- function(obj, expected_rate = 0.075, pcs = 1:20, seed = 42) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    log_msg("WARNING: DoubletFinder unavailable; skipping doublet removal. Doublets inflate")
    log_msg("         apparent NK/T 'intermediate' states -- record this in the methods.")
    obj$is_doublet <- FALSE
    return(obj)
  }
  set.seed(seed)
  sweep_fn <- if (exists("paramSweep", where = asNamespace("DoubletFinder"))) {
    DoubletFinder::paramSweep
  } else {
    get("paramSweep_v3", envir = asNamespace("DoubletFinder"))
  }
  find_fn <- if (exists("doubletFinder", where = asNamespace("DoubletFinder"))) {
    DoubletFinder::doubletFinder
  } else {
    get("doubletFinder_v3", envir = asNamespace("DoubletFinder"))
  }

  sweep_res <- sweep_fn(obj, PCs = pcs, sct = FALSE)
  stats <- DoubletFinder::summarizeSweep(sweep_res, GT = FALSE)
  pk <- as.numeric(as.character(
    DoubletFinder::find.pK(stats)$pK[which.max(DoubletFinder::find.pK(stats)$BCmetric)]
  ))
  n_expected <- round(expected_rate * ncol(obj))
  obj <- find_fn(obj, PCs = pcs, pN = 0.25, pK = pk, nExp = n_expected, sct = FALSE)
  class_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)[1]
  obj$is_doublet <- obj@meta.data[[class_col]] == "Doublet"
  log_msg("doublets: %d/%d flagged (pK=%.3f)", sum(obj$is_doublet), ncol(obj), pk)
  subset(obj, subset = is_doublet == FALSE)
}

# ------------------------------------------------------------------ integration

#' Integrate within cancer type, over patient.
#'
#' Cross-cancer-type alignment is deliberately NOT performed at the embedding
#' level. Any comparison across cancer types is made on shared cell-type labels
#' and per-sample summary statistics, so that the pan-cancer claim rests on
#' replication across cohorts rather than on a shared latent space that may have
#' erased the difference being tested.
integrate_within_cancer_type <- function(obj, patient_key = "patient",
                                         cancer_key = "cancer_type",
                                         n_pcs = 30, seed = 42) {
  stop_if(!all(c(patient_key, cancer_key) %in% colnames(obj@meta.data)),
          "metadata must contain '%s' and '%s'", patient_key, cancer_key)

  results <- list()
  for (ct in unique(obj@meta.data[[cancer_key]])) {
    log_msg("integrating cancer type %s", ct)
    sub <- subset(obj, cells = colnames(obj)[obj@meta.data[[cancer_key]] == ct])
    n_patients <- length(unique(sub@meta.data[[patient_key]]))
    if (n_patients < 2) {
      log_msg("  only %d patient(s); skipping integration for %s", n_patients, ct)
      sub <- NormalizeData(sub, verbose = FALSE) |>
        FindVariableFeatures(verbose = FALSE) |>
        ScaleData(verbose = FALSE) |>
        RunPCA(npcs = n_pcs, seed.use = seed, verbose = FALSE)
      sub[["harmony"]] <- sub[["pca"]]
      results[[ct]] <- sub
      next
    }
    sub[["RNA"]] <- split(sub[["RNA"]], f = sub@meta.data[[patient_key]])
    sub <- NormalizeData(sub, verbose = FALSE) |>
      FindVariableFeatures(verbose = FALSE) |>
      ScaleData(verbose = FALSE) |>
      RunPCA(npcs = n_pcs, seed.use = seed, verbose = FALSE)
    sub <- IntegrateLayers(sub, method = HarmonyIntegration,
                           orig.reduction = "pca", new.reduction = "harmony",
                           verbose = FALSE)
    sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
    results[[ct]] <- sub
    log_msg("  %s: %d cells, %d patients integrated", ct, ncol(sub), n_patients)
  }
  results
}

# ------------------------------------------------------------------ deconvolution

build_reference <- function(sc_obj, celltype_key = "cell_type", min_cells = 25) {
  counts <- GetAssayData(sc_obj, assay = "RNA", layer = "counts")
  labels <- factor(sc_obj@meta.data[[celltype_key]])
  keep <- labels %in% names(which(table(labels) >= min_cells))
  dropped <- setdiff(levels(labels), levels(droplevels(labels[keep])))
  if (length(dropped)) {
    log_msg("reference: dropping %s (<%d cells); they cannot be deconvolved and their",
            paste(dropped, collapse = ", "), min_cells)
    log_msg("           signal will be absorbed by the nearest retained type")
  }
  Reference(counts[, keep], droplevels(labels[keep]), colSums(counts[, keep]))
}

run_rctd <- function(st_obj, reference, mode = "full", max_cores = 4) {
  coords <- GetTissueCoordinates(st_obj)[, c("x", "y"), drop = FALSE]
  counts <- GetAssayData(st_obj, assay = "Spatial", layer = "counts")
  puck <- SpatialRNA(coords, counts, colSums(counts))
  rctd <- create.RCTD(puck, reference, max_cores = max_cores)
  rctd <- run.RCTD(rctd, doublet_mode = mode)
  weights <- normalize_weights(rctd@results$weights)
  log_msg("RCTD (%s): %d spots x %d cell types", mode, nrow(weights), ncol(weights))
  as.matrix(weights)
}

#' Benchmark deconvolution against simulated spot mixtures of known composition.
#'
#' This is the check that decides which cell types may carry a conclusion. Spots
#' are built by summing counts from a known number of reference cells, so the
#' ground-truth proportions are exact. Per-cell-type Pearson r between true and
#' estimated proportion becomes the reliability flag.
benchmark_deconvolution <- function(sc_obj, reference, celltype_key = "cell_type",
                                    n_spots = 500, cells_per_spot = 8,
                                    reliable_r = 0.7, max_cores = 4, seed = 42) {
  set.seed(seed)
  counts <- GetAssayData(sc_obj, assay = "RNA", layer = "counts")
  labels <- as.character(sc_obj@meta.data[[celltype_key]])
  types <- sort(unique(labels))

  sim <- Matrix(0, nrow = nrow(counts), ncol = n_spots, sparse = TRUE)
  truth <- matrix(0, nrow = n_spots, ncol = length(types), dimnames = list(NULL, types))
  for (i in seq_len(n_spots)) {
    picked <- sample(ncol(counts), cells_per_spot)
    sim[, i] <- Matrix::rowSums(counts[, picked, drop = FALSE])
    tab <- table(factor(labels[picked], levels = types))
    truth[i, ] <- as.numeric(tab) / cells_per_spot
  }
  rownames(sim) <- rownames(counts)
  colnames(sim) <- paste0("sim", seq_len(n_spots))

  coords <- data.frame(x = runif(n_spots), y = runif(n_spots), row.names = colnames(sim))
  puck <- SpatialRNA(coords, sim, colSums(sim))
  rctd <- run.RCTD(create.RCTD(puck, reference, max_cores = max_cores), doublet_mode = "full")
  est <- as.matrix(normalize_weights(rctd@results$weights))

  shared <- intersect(colnames(truth), colnames(est))
  report <- data.frame(
    cell_type = shared,
    true_mean_proportion = colMeans(truth[, shared, drop = FALSE]),
    pearson_r = vapply(shared, function(ct) {
      suppressWarnings(cor(truth[, ct], est[, ct]))
    }, numeric(1)),
    rmse = vapply(shared, function(ct) {
      sqrt(mean((truth[, ct] - est[, ct])^2))
    }, numeric(1)),
    row.names = NULL
  )
  report$pearson_r[is.na(report$pearson_r)] <- 0
  report$reliable <- report$pearson_r >= reliable_r

  for (i in seq_len(nrow(report))) {
    log_msg("  benchmark %-22s r=%.2f rmse=%.3f mean_prop=%.3f -> %s",
            report$cell_type[i], report$pearson_r[i], report$rmse[i],
            report$true_mean_proportion[i],
            ifelse(report$reliable[i], "reliable", "UNRELIABLE"))
  }
  unreliable <- report$cell_type[!report$reliable]
  if (length(unreliable)) {
    log_msg("WARNING: %s failed the mixture benchmark. Their weights are similarity scores,",
            paste(unreliable, collapse = ", "))
    log_msg("         not abundances, and must not carry a conclusion on this platform.")
  }
  report
}

#' Agreement between RCTD full and doublet mode, per cell type.
#'
#' Full mode always returns a weight for every type; doublet mode commits to at
#' most two. Where the two disagree for a type, the full-mode weight is mostly
#' reference similarity rather than evidence of presence.
mode_agreement <- function(weights_full, rctd_doublet_results, types) {
  called <- table(factor(c(rctd_doublet_results$first_type,
                           rctd_doublet_results$second_type), levels = types))
  data.frame(
    cell_type = types,
    full_mean_weight = colMeans(weights_full[, types, drop = FALSE]),
    doublet_call_rate = as.numeric(called) / nrow(weights_full),
    row.names = NULL
  )
}

# ------------------------------------------------------------------ entrypoint

main <- function() {
  option_list <- list(
    make_option("--sc-input", type = "character", help = "scRNA-seq Seurat object (.rds)"),
    make_option("--st-input", type = "character", help = "spatial Seurat object (.rds)"),
    make_option("--outdir", type = "character", default = "results/integration"),
    make_option("--celltype-key", type = "character", default = "cell_type"),
    make_option("--patient-key", type = "character", default = "patient"),
    make_option("--cancer-key", type = "character", default = "cancer_type"),
    make_option("--nmads", type = "double", default = 3),
    make_option("--max-cores", type = "integer", default = 4),
    make_option("--seed", type = "integer", default = 42),
    make_option("--skip-benchmark", action = "store_true", default = FALSE)
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  stop_if(is.null(opt$`sc-input`) || is.null(opt$`st-input`),
          "--sc-input and --st-input are required")
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  set.seed(opt$seed)
  check_versions()

  log_msg("loading %s", opt$`sc-input`)
  sc_obj <- readRDS(opt$`sc-input`)

  samples <- unique(sc_obj@meta.data[[opt$`patient-key`]])
  log_msg("QC over %d samples", length(samples))
  qc_log <- list()
  keep_cells <- character(0)
  for (s in samples) {
    cells <- colnames(sc_obj)[sc_obj@meta.data[[opt$`patient-key`]] == s]
    sub <- subset(sc_obj, cells = cells)
    keep <- qc_filter_sample(sub, s, nmads = opt$nmads)
    qc_log[[s]] <- attr(keep, "thresholds")
    keep_cells <- c(keep_cells, cells[keep])
  }
  sc_obj <- subset(sc_obj, cells = keep_cells)
  write.csv(do.call(rbind, lapply(qc_log, as.data.frame)),
            file.path(opt$outdir, "qc_thresholds.csv"), row.names = FALSE)

  sc_obj <- remove_doublets(sc_obj, seed = opt$seed)
  integrated <- integrate_within_cancer_type(
    sc_obj, opt$`patient-key`, opt$`cancer-key`, seed = opt$seed
  )
  saveRDS(integrated, file.path(opt$outdir, "integrated_by_cancer_type.rds"))

  log_msg("building deconvolution reference")
  reference <- build_reference(sc_obj, opt$`celltype-key`)

  if (!opt$`skip-benchmark`) {
    log_msg("benchmarking deconvolution on simulated mixtures")
    report <- benchmark_deconvolution(sc_obj, reference, opt$`celltype-key`,
                                      max_cores = opt$`max-cores`, seed = opt$seed)
    write.csv(report, file.path(opt$outdir, "deconvolution_reliability.csv"), row.names = FALSE)
  }

  log_msg("loading %s", opt$`st-input`)
  st_obj <- readRDS(opt$`st-input`)
  weights_full <- run_rctd(st_obj, reference, "full", opt$`max-cores`)
  saveRDS(weights_full, file.path(opt$outdir, "rctd_weights_full.rds"))

  log_msg("done -> %s", opt$outdir)
}

if (sys.nframe() == 0) main()
