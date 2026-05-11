# =============================================================================
# burden_analyses.R
# Cryptic splicing burden - motor cortex, frontal/temporal, RiMOD, spinal cord

# =============================================================================


# =============================================================================
# SECTION 1 - LIBRARIES & GLOBAL CONSTANTS
# =============================================================================

library(arrow)
library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(patchwork)
library(pROC)
library(rstatix)
library(pheatmap)
library(cellGeometry)
library(glue)
library(janitor)

PLOT_TEXT_SIZE  <- 14    # base ggplot text size for all plots
TOP_N_GO        <- 10    # top N GO terms per ontology in enrichment dot plots
FDR_CUTOFF      <- 0.01  # FDR threshold for defining up/down gene sets
MIN_MEAN_LOGTPM <- 0.1   # minimum mean log10-TPM to retain a gene in GLM
LOG_BASE        <- 10    # log base used throughout (log10)

log_tpm <- function(x) log(x, base = LOG_BASE)


# =============================================================================
# SECTION 2 - PATHS
# All external file references are declared here. Adjust to your setup.
# Files marked [COPY NEEDED] must be manually placed in data/ before running.
# Parquet datasets are too large for git; queried results are cached to data/.
# =============================================================================

NYGC_METADATA_PATH <- "data/metadta_nygc.csv"
# Key columns used: sample, rna_tissue_source_simplified, tdp_proteinopathy,
#                   simplified_mutations, sex_call, rin

# RiMOD metadata - two files serve different purposes:
#   rimod_meta.csv  : sample-level info incl. tdp_proteinopathy [COPY NEEDED]
#                     Copy from junction_explorer_app/data/rimod/rimod_meta.csv
#   metadata_rimod.csv : sex and RIN 
RIMOD_META_PATH         <- "data/rimod_meta.csv"
META_RIMOD_WITH_SEX_PATH <- "data/metadata_rimod.csv"

# Junction quality metrics 
SUMMARY_POST_PATH <- "data/summary_metrics_post_mortem.csv"

# Deconvolution fit objects [COPY NEEDED from Desktop]

FRONTAL_FIT_PATH <- "data/ba9_frontaltemporal_fit.RDS"
MOTOR_FIT_PATH   <- "data/ba4_motor_fit.RDS"
CORD_FIT_PATH    <- "data/cord_fit.RDS"

# Spinal cord bulk RNA

CORD_BULK_PATH     <- "data/cord_bulk_gene_counts.csv"
CORD_GENE_MAP_PATH <- "data/cord_gene_mapper.csv"

# Parquet datasets (live in junction_explorer_app)
NYGC_PARQUETS_PATH  <- "/Users/annaleigh/Documents/GitHub/junction_explorer_app/data/nygc/parquets"
RIMOD_PARQUETS_PATH <- "/Users/annaleigh/Documents/GitHub/junction_explorer_app/data/rimod/parquets"

# Queried junction caches - from above
CACHE_NYGC_PATH  <- "data/nygc_queried.rds"
CACHE_RIMOD_PATH <- "data/rimod_queried.rds"

# Output directory for saved figures
OUTPUT_DIR <- "output"
dir.create(OUTPUT_DIR, showWarnings = FALSE)


# =============================================================================
# SECTION 3 - FUNCTION DEFINITIONS
# =============================================================================

# -----------------------------------------------------------------------------
# 3a. Data loading
# -----------------------------------------------------------------------------

#' Query one junction from a parquet dataset
#' @param junc  "chr1:12345-67890" format
#' @param dataset_parquet  Path to parquet directory
query_junction <- function(junc, dataset_parquet) {
    parts <- unlist(strsplit(junc, "[:\\-]"))
    arrow::open_dataset(dataset_parquet) %>%
        filter(chrom == parts[[1]],
               start == as.numeric(parts[[2]]),
               end   == as.numeric(parts[[3]])) %>%
        collect()
}

#' Tidy a cellGeometry fit object into sub/group proportion data.frames
#' @param fit           cellGeometry fit object
#' @param sample_prefix Optional prefix to strip (e.g. "rimo" for RiMOD)
#' @return list($sub, $group) each with a 'sample' column
load_deconv <- function(fit, sample_prefix = NULL) {
    tidy_pct <- function(mat) {
        df <- mat %>% as.data.frame() %>% tibble::rownames_to_column("sample")
        if (!is.null(sample_prefix))
            df <- df %>%
                filter(grepl(sample_prefix, sample)) %>%
                tidyr::separate(sample, into = "sample", extra = "drop")
        df
    }
    list(sub   = tidy_pct(fit$subclass$percent),
         group = tidy_pct(fit$group$percent))
}


# -----------------------------------------------------------------------------
# 3b. Cryptic burden
# -----------------------------------------------------------------------------

#' Compute per-sample cryptic splicing burden
#'
#' Filters to junctions passing QC for a tissue group, z-scores PSI per
#' junction, then sums across junctions per sample.
#'
#' @param queried_data  Bound parquet rows (with chrom/start/end/psi/sample)
#' @param samples_meta  Metadata; must contain sample_col + tdp_proteinopathy
#' @param summary_post  Master junction QC table
#' @param group_label   String matching summary_post$group for this tissue
#' @param sample_col    Column in samples_meta that maps to parquet 'sample'
#' @param fpr_thresh / tpr_thresh / min_tdp_path  QC thresholds
#' @return data.frame: sample, cryptic_burden, + all samples_meta columns
compute_burden <- function(queried_data,
                           samples_meta,
                           summary_post,
                           group_label,
                           sample_col   = "sample",
                           fpr_thresh   = 0.1,
                           tpr_thresh   = 0.1,
                           min_tdp_path = 5) {

    passing <- summary_post %>%
        filter(fpr_value       <= fpr_thresh,
               tpr_value       >= tpr_thresh,
               value_tdp_path  >= min_tdp_path,
               group           == group_label) %>%
        pull(paste_into_igv_junction)

    queried_data %>%
        mutate(paste_into_igv_junction = glue("{chrom}:{start}-{end}")) %>%
        filter(sample                  %in% samples_meta[[sample_col]],
               paste_into_igv_junction %in% passing) %>%
        group_by(paste_into_igv_junction) %>%
        mutate(scale_psi = scale(psi) %>% as.vector()) %>%
        ungroup() %>%
        distinct(sample, paste_into_igv_junction, scale_psi) %>%
        group_by(sample) %>%
        summarise(cryptic_burden = sum(scale_psi), .groups = "drop") %>%
        full_join(samples_meta, by = setNames(sample_col, "sample"))
}


# -----------------------------------------------------------------------------
# 3c. ROC / pROC plots
# -----------------------------------------------------------------------------

#' Wide "top-two junction" data.frame needed by make_roc_plot
prep_top_two <- function(samples_long, top2_for_group, burden_df) {
    samples_long %>%
        left_join(burden_df %>% distinct(sample, numeric_prediction),
                  by = "sample") %>%
        filter(paste_into_igv_junction %in% top2_for_group$paste_into_igv_junction) %>%
        distinct(paste_into_igv_junction, psi, sample, numeric_prediction) %>%
        pivot_wider(names_from  = "paste_into_igv_junction",
                    values_from = "psi") %>%
        mutate(across(starts_with("chr"), ~ scale(.) %>% as.vector())) %>%
        select(where(~ !any(is.nan(.))))
}

#' Three-curve ROC plot: burden + top 2 individual junctions
#'
#' @param burden_df    Must have numeric_prediction and cryptic_burden
#' @param top_two_df   Wide df from prep_top_two()
#' @param junc1_col / junc2_col   Column names for the two junctions
#' @param gene1_label / gene2_label  Legend labels
#' @param title        Plot title
make_roc_plot <- function(burden_df,
                          top_two_df,
                          junc1_col,
                          junc2_col,
                          gene1_label,
                          gene2_label,
                          title = "") {

    roc_burden <- pROC::roc_(data = burden_df, "numeric_prediction", "cryptic_burden", quiet = TRUE)
    roc_j1     <- pROC::roc_(data = top_two_df, "numeric_prediction", junc1_col, quiet = TRUE)
    roc_j2     <- pROC::roc_(data = top_two_df, "numeric_prediction", junc2_col, quiet = TRUE)

    dev.new(width = 5, height = 5)

    plot(roc_burden, col = "#D64933", lwd = 2.3, main = title, legacy.axes = TRUE)
    plot(roc_j1,     col = "#730071", lwd = 2.3, add = TRUE)
    plot(roc_j2,     col = "#92DCE5", lwd = 2.3, add = TRUE)
    legend("bottomright",
           legend = c(
               sprintf("Cryptic Burden (AUC = %.3f)", pROC::auc(roc_burden)),
               glue::glue("{gene1_label} (AUC = {round(pROC::auc(roc_j1), 3)})"),
               glue::glue("{gene2_label} (AUC = {round(pROC::auc(roc_j2), 3)})")
           ),
           col = c("#D64933", "#730071", "#92DCE5"),
           lwd = 4, bty = "n")

    p <- recordPlot()
    dev.off()
    p
}


# -----------------------------------------------------------------------------
# 3d. Deconvolution correlation
# -----------------------------------------------------------------------------

#' Spearman correlation of cryptic burden with cell-type proportions
#'
#' Works identically for subtype and major-group tables.
#'
#' @param burden_df   Must have: sample, cryptic_burden, tdp_proteinopathy
#' @param deconv_df   Proportion df with a 'sample' column + cell-type columns
#' @param tdp_case    Filter burden_df to this disease label before correlating
#' @return data.frame with cor, p, p_signif per cell type
correlate_burden_deconv <- function(burden_df,
                                    deconv_df,
                                    tdp_case = "ALS-TDP") {
    burden_df %>%
        filter(tdp_proteinopathy == tdp_case) %>%
        left_join(deconv_df, by = "sample") %>%
        select(sample, cryptic_burden,
               where(is.numeric) & !matches("cryptic_burden|numeric_prediction")) %>%
        reshape2::melt(id.vars = c("sample", "cryptic_burden")) %>%
        group_by(variable) %>%
        rstatix::cor_test(cryptic_burden, value, method = "spearman") %>%
        mutate(p_signif = case_when(
            p < 0.001 ~ "***",
            p < 0.01  ~ "**",
            p < 0.05  ~ "*",
            TRUE      ~ ""
        ))
}

#' Bar chart of deconvolution-burden correlations
#'
#' @param cor_df       Output of correlate_burden_deconv()
#' @param cell_table   Named vector: cell name -> class (from fit$mk$cell_table).
#'                     Pass NULL to skip class-level colouring.
#' @param title        Plot title
plot_deconv_cor <- function(cor_df,
                            cell_table    = NULL,
                            class_colours = c(In   = "#2E6F95",
                                              Ex   = "#E64A19",
                                              Glia = "#9B6BBA",
                                              Vasc = "#A32C34"),
                            title = "") {
    if (!is.null(cell_table)) {
        cor_df <- cor_df %>%
            left_join(tibble::enframe(cell_table,
                                      name  = "variable",
                                      value = "cell_class"),
                      by = "variable") %>%
            mutate(cell_class = as.character(cell_class)) %>%
            mutate(cell_class = ifelse(is.na(cell_class), variable, cell_class))

        fill_aes   <- aes(fill = cell_class)
        fill_scale <- scale_fill_manual(values = class_colours, na.value = "grey70")
    } else {
        cor_df$cell_class <- "unknown"
        fill_aes   <- aes(fill = cell_class)
        fill_scale <- scale_fill_manual(values = c(unknown = "steelblue"))
    }

    ggplot(cor_df, aes(x = reorder(variable, cor), y = cor)) +
        geom_col(fill_aes) +
        geom_text(aes(label = p_signif), size = 5, color = "black") +
        coord_flip() +
        fill_scale +
        theme_bw(base_size = PLOT_TEXT_SIZE) +
        labs(x = "", y = "Spearman ρ with cryptic burden",
             fill = "Cell class", title = title) +
        theme(panel.grid.major.y = element_blank())
}


# =============================================================================
# SECTION 4 - LOAD SHARED DATA
# =============================================================================

# --- Metadata ----------------------------------------------------------------
nygc_metadata       <- fread(NYGC_METADATA_PATH)
rimod_meta          <- fread(RIMOD_META_PATH)
meta_rimod_with_sex <- fread(META_RIMOD_WITH_SEX_PATH)
summary_post        <- fread(SUMMARY_POST_PATH) %>%
    mutate(junctions_coords = paste_into_igv_junction)

# --- Query passing junctions (cached after first run) ------------------------
ny_junctions <- summary_post %>%
    filter(fpr_value <= 0.1, tpr_value >= 0.1,
           value_tdp_path >= 5, group != "RiMOD") %>%
    distinct(paste_into_igv_junction)

ri_junctions <- summary_post %>%
    filter(fpr_value <= 0.1, tpr_value >= 0.1,
           value_tdp_path >= 5, group == "RiMOD") %>%
    distinct(paste_into_igv_junction)

if (file.exists(CACHE_NYGC_PATH)) {
    message("Loading cached NYGC junction data from ", CACHE_NYGC_PATH)
    nygc_queried <- readRDS(CACHE_NYGC_PATH)
} else {
    message("Querying NYGC parquets (this may take a while)...")
    nygc_queried <- purrr::map(ny_junctions$paste_into_igv_junction,
                               query_junction, NYGC_PARQUETS_PATH) %>%
        rbindlist() %>%
        mutate(paste_into_igv_junction = glue("{chrom}:{start}-{end}"))
    saveRDS(nygc_queried, CACHE_NYGC_PATH)
    message("Cached NYGC query to ", CACHE_NYGC_PATH)
}

if (file.exists(CACHE_RIMOD_PATH)) {
    message("Loading cached RiMOD junction data from ", CACHE_RIMOD_PATH)
    rimod_queried <- readRDS(CACHE_RIMOD_PATH)
} else {
    message("Querying RiMOD parquets (this may take a while)...")
    rimod_queried <- purrr::map(ri_junctions$paste_into_igv_junction,
                                query_junction, RIMOD_PARQUETS_PATH) %>%
        rbindlist() %>%
        mutate(paste_into_igv_junction = glue("{chrom}:{start}-{end}"),
               psi = junction_count / cluster_count,
               psi = replace_na(psi, 0))
    saveRDS(rimod_queried, CACHE_RIMOD_PATH)
    message("Cached RiMOD query to ", CACHE_RIMOD_PATH)
}

# Top-2 junctions per tissue (for ROC plots)
top2 <- summary_post %>%
    filter(fpr_value <= 0.1, tpr_value >= 0.1, value_tdp_path >= 5) %>%
    group_by(group) %>%
    slice_max(auc, n = 2) %>%
    as.data.table()


# =============================================================================
# SECTION 5 - DECONVOLUTION FITS
# =============================================================================

frontal_fit    <- readRDS(FRONTAL_FIT_PATH)
motor_fit      <- readRDS(MOTOR_FIT_PATH)
spinal_fit_raw <- readRDS(CORD_FIT_PATH)

spinal_mk <- updateMarkers(
    spinal_fit_raw$mk,
    remove_subclass = c("Microglia-3", "Oligo-5", "Oligo-2",
                        "Oligo-6", "Ex-Dorsal-9", "Ex-V-2")
)

cord_gene_mapper <- read.csv(CORD_GENE_MAP_PATH)
cord_bulk <- fread(CORD_BULK_PATH) %>%
    left_join(cord_gene_mapper %>% distinct(cord_name, tpm_calc_name),
              by = c("V1" = "tpm_calc_name")) %>%
    filter(!is.na(cord_name)) %>%
    tibble::column_to_rownames("cord_name") %>%
    select(-V1)

spinal_fit <- deconvolute(spinal_mk, cord_bulk)

# Tidy deconvolution proportions
deconv <- list(
    motor   = load_deconv(motor_fit),
    frontal = load_deconv(frontal_fit),
    rimod   = load_deconv(frontal_fit, sample_prefix = "rimo"),
    cord    = load_deconv(spinal_fit)
)


# =============================================================================
# SECTION 6 - SAMPLE METADATA SUBSETS
# =============================================================================

# NYGC metadata uses 'sample' as the parquet sample identifier column.
motor_meta <- nygc_metadata %>%
    filter(rna_tissue_source_simplified == "Cortex_Motor",
           tdp_proteinopathy %in% c("Control", "ALS-non-TDP", "ALS-TDP")) %>%
    distinct(rna_tissue_source_simplified, sample, tdp_proteinopathy, simplified_mutations)

frontal_meta <- nygc_metadata %>%
    filter(rna_tissue_source_simplified %in% c("Cortex_Frontal", "Cortex_Temporal"),
           tdp_proteinopathy %in% c("Control", "FTLD-non-TDP", "FTLD-TDP")) %>%
    distinct(rna_tissue_source_simplified, sample, tdp_proteinopathy, simplified_mutations)

cord_meta <- nygc_metadata %>%
    filter(grepl("Cord", rna_tissue_source_simplified),
           tdp_proteinopathy %in% c("Control", "ALS-non-TDP", "ALS-TDP")) %>%
    distinct(rna_tissue_source_simplified, sample, tdp_proteinopathy, simplified_mutations)

# rimod_meta loaded above; sex joined in pca_go_analyses.R GLM section


# =============================================================================
# SECTION 7 - CRYPTIC BURDEN (ALL DATASETS)
# =============================================================================

motor_burden   <- compute_burden(nygc_queried,  motor_meta,   summary_post, "Motor Cortex",        sample_col = "sample")
frontal_burden <- compute_burden(nygc_queried,  frontal_meta, summary_post, "FTD frontal/temporal", sample_col = "sample")
rimod_burden   <- compute_burden(rimod_queried, rimod_meta,   summary_post, "RiMOD",               sample_col = "sample")
cord_burden    <- compute_burden(nygc_queried,  cord_meta,    summary_post, "Cord",                sample_col = "sample")

# Add binary outcome column for ROC
motor_burden   <- motor_burden   %>% mutate(numeric_prediction = as.integer(tdp_proteinopathy == "ALS-TDP"))
frontal_burden <- frontal_burden %>% mutate(numeric_prediction = as.integer(tdp_proteinopathy == "FTLD-TDP"))
rimod_burden   <- rimod_burden   %>% mutate(numeric_prediction = as.integer(tdp_proteinopathy == "FTLD-TDP"))
cord_burden    <- cord_burden    %>% mutate(numeric_prediction = as.integer(tdp_proteinopathy == "ALS-TDP"))


# =============================================================================
# SECTION 8 - LONG-FORMAT DATA (required by ROC and pca_go_analyses.R)
# =============================================================================

#' Join queried junction data back to metadata with a common 'sample' column
join_meta <- function(queried, meta, sample_col) {
    these_samples <- meta %>% pull(!!sample_col)
    meta_renamed  <- meta %>% dplyr::rename(sample = !!sample_col)
    queried %>%
        filter(sample %in% these_samples) %>%
        left_join(meta_renamed, by = "sample")
}

motor_long   <- join_meta(nygc_queried,  motor_meta,   "sample")
frontal_long <- join_meta(nygc_queried,  frontal_meta, "sample")
rimod_long   <- join_meta(rimod_queried, rimod_meta,   "sample")
cord_long    <- join_meta(nygc_queried,  cord_meta,    "sample")


# =============================================================================
# SECTION 9 - ROC PLOTS (ALL DATASETS)
# =============================================================================

# Motor cortex ----------------------------------------------------------------
motor_top2    <- top2 %>% filter(group == "Motor Cortex")
motor_top_two <- prep_top_two(motor_long, motor_top2, motor_burden)
roc_motor     <- make_roc_plot(
    motor_burden, motor_top_two,
    junc1_col   = motor_top2 %>% slice_max(auc) %>% pull(paste_into_igv_junction),
    junc2_col   = motor_top2 %>% slice_min(auc) %>% pull(paste_into_igv_junction),
    gene1_label = motor_top2 %>% slice_max(auc) %>% pull(gene_name),
    gene2_label = motor_top2 %>% slice_min(auc) %>% pull(gene_name),
    title       = "NYGC Motor Cortex"
)

# Frontal/temporal ------------------------------------------------------------
frontal_top2    <- top2 %>% filter(group == "FTD frontal/temporal")
frontal_top_two <- prep_top_two(frontal_long, frontal_top2, frontal_burden)
roc_frontal     <- make_roc_plot(
    frontal_burden, frontal_top_two,
    junc1_col   = frontal_top2 %>% slice_max(auc) %>% pull(paste_into_igv_junction),
    junc2_col   = frontal_top2 %>% slice_min(auc) %>% pull(paste_into_igv_junction),
    gene1_label = frontal_top2 %>% slice_max(auc) %>% pull(gene_name),
    gene2_label = frontal_top2 %>% slice_min(auc) %>% pull(gene_name),
    title       = "NYGC Frontal/Temporal"
)

# RiMOD -----------------------------------------------------------------------
rimod_top2    <- top2 %>% filter(group == "RiMOD")
rimod_top_two <- prep_top_two(rimod_long, rimod_top2, rimod_burden)
roc_rimod     <- make_roc_plot(
    rimod_burden, rimod_top_two,
    junc1_col   = rimod_top2 %>% slice_max(auc) %>% pull(paste_into_igv_junction),
    junc2_col   = rimod_top2 %>% slice_min(auc) %>% pull(paste_into_igv_junction),
    gene1_label = rimod_top2 %>% slice_max(auc) %>% pull(gene_name),
    gene2_label = rimod_top2 %>% slice_min(auc) %>% pull(gene_name),
    title       = "RiMOD Frontal Cortex"
)

# Spinal cord -----------------------------------------------------------------
cord_top2    <- top2 %>% filter(group == "Cord")
cord_top_two <- prep_top_two(cord_long, cord_top2, cord_burden)
roc_cord     <- make_roc_plot(
    cord_burden, cord_top_two,
    junc1_col   = cord_top2 %>% slice_max(auc) %>% pull(paste_into_igv_junction),
    junc2_col   = cord_top2 %>% slice_min(auc) %>% pull(paste_into_igv_junction),
    gene1_label = cord_top2 %>% slice_max(auc) %>% pull(gene_name),
    gene2_label = cord_top2 %>% slice_min(auc) %>% pull(gene_name),
    title       = "NYGC Spinal Cord"
)


# =============================================================================
# SECTION 10 - DECONVOLUTION CORRELATIONS (ALL DATASETS)
# =============================================================================

# Motor cortex ----------------------------------------------------------------
motor_cor_sub   <- correlate_burden_deconv(motor_burden,   deconv$motor$sub,   tdp_case = "ALS-TDP")
motor_cor_major <- correlate_burden_deconv(motor_burden,   deconv$motor$group, tdp_case = "ALS-TDP")

plot_deconv_cor(motor_cor_sub,
                cell_table = frontal_fit$mk$cell_table,
                title      = "ALS-TDP motor cortex - subtype deconvolution")
plot_deconv_cor(motor_cor_major,
                cell_table = frontal_fit$mk$cell_table,
                title      = "ALS-TDP motor cortex - major cell groups")

# Frontal/temporal ------------------------------------------------------------
frontal_cor_sub   <- correlate_burden_deconv(frontal_burden, deconv$frontal$sub,   tdp_case = "FTLD-TDP")
frontal_cor_major <- correlate_burden_deconv(frontal_burden, deconv$frontal$group, tdp_case = "FTLD-TDP")

plot_deconv_cor(frontal_cor_sub,
                cell_table = frontal_fit$mk$cell_table,
                title      = "FTLD-TDP frontal/temporal - subtype deconvolution")
plot_deconv_cor(frontal_cor_major,
                cell_table = frontal_fit$mk$cell_table,
                title      = "FTLD-TDP frontal/temporal - major cell groups")

# RiMOD -----------------------------------------------------------------------
rimod_cor_sub   <- correlate_burden_deconv(
    rimod_burden %>% select(sample, cryptic_burden, tdp_proteinopathy),
    deconv$rimod$sub,   tdp_case = "FTLD-TDP")
rimod_cor_major <- correlate_burden_deconv(
    rimod_burden %>% select(sample, cryptic_burden, tdp_proteinopathy),
    deconv$rimod$group, tdp_case = "FTLD-TDP")

plot_deconv_cor(rimod_cor_sub,
                cell_table = frontal_fit$mk$cell_table,
                title      = "FTLD-TDP RiMOD - subtype deconvolution")
plot_deconv_cor(rimod_cor_major,
                cell_table = frontal_fit$mk$cell_table,
                title      = "FTLD-TDP RiMOD - major cell groups")

# Spinal cord -----------------------------------------------------------------
cord_cor_sub   <- correlate_burden_deconv(cord_burden, deconv$cord$sub,   tdp_case = "ALS-TDP")
cord_cor_major <- correlate_burden_deconv(cord_burden, deconv$cord$group, tdp_case = "ALS-TDP")

plot_deconv_cor(cord_cor_sub,
                cell_table = spinal_fit$mk$cell_table,
                title      = "ALS-TDP spinal cord - subtype deconvolution")
plot_deconv_cor(cord_cor_major,
                cell_table = spinal_fit$mk$cell_table,
                title      = "ALS-TDP spinal cord - major cell groups")


# =============================================================================
# SECTION 11 - COMBINED FIGURES
# =============================================================================


# --- Cortical deconvolution heatmap ------------------------------------------
cortical_cor_df <- bind_rows(
    motor_cor_sub   %>% mutate(group = "NYGC Motor Cortex"),
    frontal_cor_sub %>% mutate(group = "NYGC Frontal/Temporal"),
    rimod_cor_sub   %>% mutate(group = "RiMOD Frontal")
)

cor_mat <- cortical_cor_df %>%
    select(variable, group, cor) %>%
    pivot_wider(names_from = group, values_from = cor) %>%
    column_to_rownames("variable") %>%
    as.matrix()

stars_mat <- cortical_cor_df %>%
    select(variable, group, p_signif) %>%
    pivot_wider(names_from = group, values_from = p_signif) %>%
    column_to_rownames("variable") %>%
    as.matrix()

annot_row <- frontal_fit$mk$cell_table %>%
    tibble::enframe(name = "cell", value = "Cell group") %>%
    tibble::column_to_rownames("cell")

cortical_correlation_heatmap = pheatmap(
    cor_mat,
    color             = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_rows      = TRUE, cluster_cols = TRUE,
    annotation_row    = annot_row,
    annotation_colors = list(`Cell group` = c(In   = "#2E6F95", Ex   = "#E64A19",
                                              Glia = "#9B6BBA", Vasc = "#A32C34")),
    na_col            = "grey90",
    display_numbers   = stars_mat,
    number_color      = "black", fontsize_number = 16,
    main              = "Cryptic burden vs cell-type deconvolution (cortical)"
)

# --- Spinal cord deconvolution heatmap ---------------------------------------
cord_cor_mat <- cord_cor_sub %>%
    mutate(group = "NYGC Spinal Cord") %>%
    select(variable, group, cor) %>%
    pivot_wider(names_from = group, values_from = cor) %>%
    column_to_rownames("variable") %>%
    as.matrix()

cord_stars_mat <- cord_cor_sub %>%
    mutate(group = "NYGC Spinal Cord") %>%
    select(variable, group, p_signif) %>%
    pivot_wider(names_from = group, values_from = p_signif) %>%
    column_to_rownames("variable") %>%
    as.matrix()

cord_annot_row <- spinal_fit$mk$cell_table %>%
    tibble::enframe(name = "cell", value = "Cell group") %>%
    tibble::column_to_rownames("cell")

cord_correlation_heatmap = pheatmap(
    cord_cor_mat,
    color             = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_rows      = TRUE, cluster_cols = FALSE,   # single column - don't cluster cols
    annotation_row    = cord_annot_row,
    # annotation_colors = list(`Cell group` = c(In   = "#2E6F95", Ex   = "#E64A19",
    #                                           Glia = "#9B6BBA", Vasc = "#A32C34")),
    na_col            = "grey90",
    display_numbers   = cord_stars_mat,
    number_color      = "black", fontsize_number = 16,
    main              = "Cryptic burden vs cell-type deconvolution (spinal cord)"
)



