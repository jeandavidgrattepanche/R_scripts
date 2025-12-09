# ============================================================
# pipeline_microbiome.R
# PCoA (2x2) + SpiecEasi network across Bac / Pico / Nano / Micro
# ============================================================

# ---- 0) User params ---------------------------------------
# Project root (adjust if needed)
setwd("/Volumes/JayDiii2SSD/Bioinfo_2025/R_2025")

# Input tables (from our previous runs)
BAC_OTU   <- "Data/B303ing_bac_final_OTUtable.tsv"
BAC_ENV   <- "Data/B303_env_bacting.txt"
EUK_OTU   <- "Data/B303ing_euk_final_OTUtable.tsv"
EUK_ENV   <- "Data/B303_env_euking.txt"

# Column selections (as before)
BAC_COLS  <- list(b0 = 0:12, b1 = 12:31, b2 = 33:47)   # bacteria: cbind two column ranges
EUK_COLS  <- 12:122                         # eukaryotes (all ingestion cols)

# Output dir
OUTDIR    <- "plots"
dir.create(OUTDIR, showWarnings = FALSE)

# PCoA aesthetics
COLOR_AESTH <- "latitude"
SHAPE_AESTH <- "Cruise"

# Network filtering thresholds (tune as needed)
MIN_PREV <- 0.10   # keep taxa present in >= 10% samples
MIN_TOT  <- 50     # and total counts >= 50

# Labeling for graphs
TOP_K_LABELS <- 40     # how many node labels to repel/keep
SEED         <- 1

# ============================================================
# ---- 1) Libraries ------------------------------------------
suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(phyloseq)
  library(phyloseqCompanion)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(plyr)
  library(SpiecEasi)
  library(igraph)
  library(ggraph)
  library(ggrepel)
})

set.seed(SEED)

# ============================================================
# ---- 2) Helpers --------------------------------------------

.clean_ps <- function(ps) {
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  ps
}

.var_explained <- function(ord) {
  eig <- ord$values$Eigenvalues
  round(100 * eig / sum(eig), 1)
}

.make_pcoa_plot <- function(ps, ord, title) {
  ve <- .var_explained(ord)
  plot_ordination(ps, ord, color = COLOR_AESTH, shape = SHAPE_AESTH) +
    geom_point(size = 3) +
    labs(
      title = title,
      subtitle = paste0("PC1: ", ve[1], "%,  PC2: ", ve[2], "% variance explained"),
      x = paste0("PC1 (", ve[1], "%)"),
      y = paste0("PC2 (", ve[2], "%)")
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, face = "italic"))
}

.pcoa_limits <- function(ord, pad = 0.05) {
  xy <- ord$vectors[, 1:2, drop = FALSE]
  rx <- range(xy[, 1], na.rm = TRUE); ry <- range(xy[, 2], na.rm = TRUE)
  dx <- diff(rx); dy <- diff(ry)
  c(xmin = rx[1] - pad*dx, xmax = rx[2] + pad*dx,
    ymin = ry[1] - pad*dy, ymax = ry[2] + pad*dy)
}

.fix_axes <- function(p, lims) {
  p + coord_equal() +
    xlim(lims["xmin"], lims["xmax"]) +
    ylim(lims["ymin"], lims["ymax"])
}

prefix_taxa <- function(ps, prefix) {
  taxa_names(ps) <- paste0(prefix, taxa_names(ps))
  ps
}

# relabel samples to station_depth
relabel_to_station <- function(ps, station_col = "wsamples") {
  st <- as.character(sample_data(ps)[[station_col]])
  # create a unique name per sample while starting with Station for alignment
  new_names <- make.unique(st, sep = "_rep")
  sample_names(ps) <- new_names
  # keep the Station itself as a separate column for later joins
  ps
}


# Filter by prevalence + total counts
filter_ps <- function(ps, min_prev = MIN_PREV, min_tot = MIN_TOT) {
  ta <- otu_table(ps)
  by_tax <- if (taxa_are_rows(ps)) 1 else 2
  prev <- apply(ta > 0, by_tax, mean)
  tot  <- apply(ta,     by_tax, sum)
  keep <- names(prev[prev >= min_prev & tot >= min_tot])
  prune_taxa(keep, ps)
}

# Build layout + label list
pick_labels <- function(g, k = TOP_K_LABELS, metric = c("degree","betweenness")) {
  metric <- match.arg(metric)
  score <- switch(metric,
                  degree = degree(g),
                  betweenness = betweenness(g, normalized = TRUE))
  names(sort(score, decreasing = TRUE))[seq_len(min(k, vcount(g)))]
}

# Nice group palette
group_cols <- c(Bac="#1b9e77", Pico="#d95f02", Nano="#7570b3", Micro="#66a61e")

# ============================================================
# ---- 3) Load data & build phyloseq objects -----------------

# Bacteria
dat_bac <- read.table(BAC_OTU, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
otu_bac <- otu_table(cbind(dat_bac[, BAC_COLS$b1], dat_bac[, BAC_COLS$b2]),
                     taxa_are_rows = TRUE, errorIfNULL = TRUE)
env_bac <- read.table(BAC_ENV, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax_bac <- as.matrix(dat_bac[ , BAC_COLS$b0])  # Potential tax_table input
# testb <- as.matrix(data.all[ , 0:15])  # Potential tax_table input
# TAX <- tax_table(testb)
ps_bac  <- .clean_ps(phyloseq(otu_bac, sample_data(env_bac),tax_table(tax_bac)))

# Eukaryotes (then split)
dat_euk <- read.table(EUK_OTU, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
otu_euk <- otu_table(dat_euk[, EUK_COLS], taxa_are_rows = TRUE, errorIfNULL = TRUE)
env_euk <- read.table(EUK_ENV, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax_euk <- as.matrix(dat_euk[ , BAC_COLS$b0])  # Potential tax_table input
ps_euk  <- .clean_ps(phyloseq(otu_euk, sample_data(env_euk),tax_table(tax_euk)))

ps_pico  <- .clean_ps(subset_samples(ps_euk, resized == "picoplankton"))
ps_nano  <- .clean_ps(subset_samples(ps_euk, resized == "nanoplankton"))
ps_micro <- .clean_ps(subset_samples(ps_euk, resized == "microplankton"))

ps_bac_r   <- relabel_to_station(ps_bac, "wsamples")
ps_pico_r  <- relabel_to_station(ps_pico, "wsamples")
ps_nano_r  <- relabel_to_station(ps_nano, "wsamples")
ps_micro_r <- relabel_to_station(ps_micro, "wsamples")


# ============================================================
# ---- 4) PCoA 2×2 figure -----------------------------------

# Ordinations
ord_bac   <- ordinate(ps_bac_r,   method = "PCoA", distance = "bray")
ord_pico  <- ordinate(ps_pico_r,  method = "PCoA", distance = "bray")
ord_nano  <- ordinate(ps_nano_r,  method = "PCoA", distance = "bray")
ord_micro <- ordinate(ps_micro_r, method = "PCoA", distance = "bray")

# Plots
p_bac   <- .make_pcoa_plot(ps_bac_r,   ord_bac,   "Bacteria")
p_pico  <- .make_pcoa_plot(ps_pico_r,  ord_pico,  "Picoeukaryotes")
p_nano  <- .make_pcoa_plot(ps_nano_r,  ord_nano,  "Nanoeukaryotes")
p_micro <- .make_pcoa_plot(ps_micro_r, ord_micro, "Microeukaryotes")

# Global axis limits
lims_mat <- rbind(.pcoa_limits(ord_bac),
                  .pcoa_limits(ord_pico),
                  .pcoa_limits(ord_nano),
                  .pcoa_limits(ord_micro))
global_lims <- c(
  xmin = min(lims_mat[, "xmin"]), xmax = max(lims_mat[, "xmax"]),
  ymin = min(lims_mat[, "ymin"]), ymax = max(lims_mat[, "ymax"])
)

p_bac   <- .fix_axes(p_bac,   global_lims)
p_pico  <- .fix_axes(p_pico,  global_lims)
p_nano  <- .fix_axes(p_nano,  global_lims)
p_micro <- .fix_axes(p_micro, global_lims)

combined <- (p_bac | p_pico) / (p_nano | p_micro) +
  plot_annotation(
    title = "PCoA (Bray–Curtis) - Bacteria & Eukaryote Size Fractions",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  ) &
  theme(legend.position = "bottom")

combined <- combined + plot_layout(guides = "collect")

ggsave(file.path(OUTDIR, "B303ing_PCoA_4panel_sharedLegend.pdf"),
       combined, width = 12, height = 10, units = "in", device = cairo_pdf)
ggsave(file.path(OUTDIR, "B303ing_PCoA_4panel_sharedLegend.png"),
       combined, width = 12, height = 10, units = "in", dpi = 300)

# ============================================================
# ---- 5) SpiecEasi network (full + cross-group) -------------



# Intersect samples across all groups
s_int <- Reduce(intersect, list(sample_data(ps_bac_r)$wsamples,
                                sample_data(ps_pico_r)$wsamples,
                                sample_data(ps_nano_r)$wsamples,
                                sample_data(ps_micro_r)$wsamples))
ps_bac_i   <- subset_samples(ps_bac_r, wsamples %in% s_int)
ps_pico_i  <- subset_samples(ps_pico_r, wsamples %in% s_int)
ps_nano_i  <- subset_samples(ps_nano_r, wsamples %in% s_int)
ps_micro_i <- subset_samples(ps_micro_r, wsamples %in% s_int)

# Prefix taxa to track origin
ps_bac_i   <- prefix_taxa(ps_bac_i,   "Bac_")
ps_pico_i  <- prefix_taxa(ps_pico_i,  "Pico_")
ps_nano_i  <- prefix_taxa(ps_nano_i,  "Nano_")
ps_micro_i <- prefix_taxa(ps_micro_i, "Micro_")

# Merge + initial cleanup
ps_all <- merge_phyloseq(ps_bac_i, ps_pico_i, ps_nano_i, ps_micro_i)
ps_all <- .clean_ps(ps_all)

# Filter by prevalence/abundance
ps_all_f <- filter_ps(ps_all, MIN_PREV, MIN_TOT)

# Fit SpiecEasi (MB method; tune as needed)
se <- spiec.easi(ps_all_f, method = "mb",
                 lambda.min.ratio = 1e-2, nlambda = 20,
                 pulsar.params = list(rep.num = 50, seed = SEED),
                 verbose = TRUE)

adj <- getRefit(se)
# Ensure names on adjacency
if (is.null(rownames(adj))) {
  rn <- taxa_names(ps_all_f)
  rownames(adj) <- colnames(adj) <- rn
}

ig <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
if (is.null(V(ig)$name)) V(ig)$name <- rownames(adj)

# Add attributes: Group, OTU (unprefixed), Taxonomy (if available)
V(ig)$Group <- sub("^([^_]+)_.*", "\\1", V(ig)$name)
V(ig)$OTU   <- sub("^[^_]+_", "", V(ig)$name)

if (!is.null(tax_table(ps_all_f, errorIfNULL = FALSE))) {
  tx <- as.data.frame(tax_table(ps_all_f))
  tx$TaxLabel <- apply(tx, 1, function(x) paste(na.omit(x), collapse = ";"))
  # Match by unprefixed OTU ID
  m <- match(V(ig)$OTU, rownames(tx))
  V(ig)$Taxonomy <- tx$TaxLabel[m]
}

# Mark cross-group edges => play with edge to ajust amount of data on network
Vmap <- setNames(V(ig)$Group, V(ig)$name)
E(ig)$CrossGroup <- apply(ends(ig, E(ig)), 1, function(p) Vmap[p[1]] != Vmap[p[2]])

# Save full graphs/tables
write_graph(ig, file.path(OUTDIR, "spieceasi_full.graphml"), format = "graphml")
write.csv(as_data_frame(ig, what = "edges"),    file.path(OUTDIR, "spieceasi_full_edges.csv"),    row.names = FALSE)
write.csv(as_data_frame(ig, what = "vertices"), file.path(OUTDIR, "spieceasi_full_nodes.csv"),    row.names = FALSE)

# Cross-group subgraph + labeled plot with repel
g_cross <- subgraph_from_edges(ig, E(ig)[CrossGroup])

lay <- create_layout(g_cross, layout = "fr")   # attributes carry into layout
labs <- pick_labels(g_cross, k = TOP_K_LABELS, metric = "degree")

p_net <- ggraph(lay) +
  geom_edge_link(alpha = 0.3, colour = "grey65") +
  geom_node_point(aes(color = Group), size = 2.5) +
  ggrepel::geom_text_repel(
    data = subset(lay, name %in% labs),
    aes(x = x, y = y, label = OTU),
    size = 3, max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2
  ) +
  scale_color_manual(values = group_cols, name = "Origin") +
  theme_void() +
  ggtitle(sprintf("Cross-group network — labels on top %d nodes", length(labs)))

ggsave(file.path(OUTDIR, "spieceasi_cross_group.png"), p_net, width = 8, height = 7, dpi = 300)
ggsave(file.path(OUTDIR, "spieceasi_cross_group.pdf"), p_net, width = 8, height = 7, device = cairo_pdf)

# ============================================================
# ---- 6) Leave-some-origin-out (LOO) refit -------------------

loo_fit <- function(ps_all, exclude = c("Bac","Pico","Nano","Micro"),
                    min_prev = MIN_PREV, min_tot = MIN_TOT,
                    method = "mb", nlambda = 20, seed = SEED) {
 # allow vector input
  if (!all(exclude %in% c("Bac","Pico","Nano","Micro"))) {
    stop("exclude must be subset of {Bac, Pico, Nano, Micro}")
  }

  # Build regex prefix pattern like "^(Micro_|Pico_)"
  patt <- paste0("^(", paste0(exclude, collapse="|"), ")_")
  keep_taxa <- taxa_names(ps_all)[!grepl(patt, taxa_names(ps_all))]
  ps <- prune_taxa(keep_taxa, ps_all)
  ps <- filter_ps(ps, min_prev, min_tot)
  set.seed(seed)
  se <- spiec.easi(ps, method = method, lambda.min.ratio = 1e-2, nlambda = nlambda,
                   pulsar.params = list(rep.num = 50, seed = seed))
  adj <- getRefit(se)
  if (is.null(rownames(adj))) {
    rn <- taxa_names(ps); rownames(adj) <- colnames(adj) <- rn
  }
  gi <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  if (is.null(V(gi)$name)) V(gi)$name <- rownames(adj)
  V(gi)$Group <- sub("^([^_]+)_.*", "\\1", V(gi)$name)
  V(gi)$OTU   <- sub("^[^_]+_", "", V(gi)$name)
  Vmap <- setNames(V(gi)$Group, V(gi)$name)
  E(gi)$CrossGroup <- apply(ends(gi, E(gi)), 1, function(p) Vmap[p[1]] != Vmap[p[2]])
  gi
}

# Example: refit excluding Micro and Pico
g_skip_bac <- loo_fit(ps_all, exclude = c("Micro","Pico"))
write_graph(g_skip_bac, file.path(OUTDIR, "spieceasi_skip_Micro.graphml"), format = "graphml")

lay <- create_layout(g_skip_bac, layout = "fr")   # attributes carry into layout
labs <- pick_labels(g_skip_bac, k = TOP_K_LABELS, metric = "degree")

p_net1 <- ggraph(lay) +
  geom_edge_link(alpha = 0.3, colour = "grey65") +
  geom_node_point(aes(color = Group), size = 2.5) +
  ggrepel::geom_text_repel(
    data = subset(lay, name %in% labs),
    aes(x = x, y = y, label = OTU),
    size = 3, max.overlaps = Inf, box.padding = 0.3, point.padding = 0.2
  ) +
  scale_color_manual(values = group_cols, name = "Origin") +
  theme_void() +
  ggtitle(sprintf("Cross-group network - labels on top %d nodes - No MicroPico", length(labs)))

# Save as PDF
ggsave(file.path(OUTDIR, "spieceasi_skip_Micro_pico.pdf"), plot = p_net1, width = 8, height = 7)

# Save as PNG (high res)
ggsave(file.path(OUTDIR, "spieceasi_skip_Micro_pico.png"), plot = p_net1, width = 8, height = 7, dpi = 300)

# ============================================================
# ---- 7) (Optional) Spearman + FDR fallback -----------------
# Use if SpiecEasi is unavailable; works on filtered ps_all_f

# Uncomment to run:
# library(Hmisc)
# X <- t(as.matrix(otu_table(ps_all_f)))  # samples x taxa (counts; CLR not required for Spearman, but OK to use)
# rc <- Hmisc::rcorr(as.matrix(X), type = "spearman")
# R  <- rc$r; P <- rc$P
# diag(R) <- 0; diag(P) <- 1
# iu <- which(upper.tri(R), arr.ind = TRUE)
# edges <- data.frame(
#   from = colnames(R)[iu[,1]],
#   to   = colnames(R)[iu[,2]],
#   r    = R[iu],
#   p    = P[iu]
# )
# edges$q <- p.adjust(edges$p, "BH")
# edges <- subset(edges, abs(r) >= 0.5 & q <= 0.05)
# g_spear <- graph_from_data_frame(edges, directed = FALSE)
# V(g_spear)$Group <- sub("^([^_]+)_.*", "\\1", V(g_spear)$name)
# V(g_spear)$OTU   <- sub("^[^_]+_", "", V(g_spear)$name)
# write_graph(g_spear, file.path(OUTDIR, "spearman_network.graphml"), format = "graphml")

# ============================================================
message("Pipeline complete. Outputs written to: ", normalizePath(OUTDIR))
