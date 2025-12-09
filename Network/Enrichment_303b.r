# enrichment v3 sample-level

setwd("/Volumes/JayDiii2SSD/Bioinfo_2025/R_2025")
# Output dir
OUTDIR    <- "plots/Nov25"

## Prepare OTU matrices and metadata
# helper to get OTU matrix with taxa as rows
get_otu_mat <- function(ps) {
  m <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) m <- t(m)
  m
}

X_bac <- get_otu_mat(ps_bac_filt)
X_euk <- get_otu_mat(ps_euk_filt)

# ensure we have a SampleID column in meta2
if (!"SampleID" %in% colnames(meta2)) {
  meta2$SampleID <- rownames(meta2)
}

# keep only samples present in both matrices and metadata
common_sids <- Reduce(
  intersect,
  list(colnames(X_bac), colnames(X_euk), meta2$SampleID)
)

X_bac <- X_bac[, common_sids, drop = FALSE]
X_euk <- X_euk[, common_sids, drop = FALSE]
meta2 <- meta2[match(common_sids, meta2$SampleID), ]


## Prepare OTU matrices and metadata
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

# unique significant pairs with their Season / NorthSouth stratum
edges_pairs <- edges_SR_annot %>%
  dplyr::mutate(pair = paste(Bac, Euk, sep = "||")) %>%
  dplyr::select(pair, Bac, Euk, Season, NorthSouth) %>%
  dplyr::distinct()

# function to get sample-level presence for one edge row
get_pair_sample_presence <- function(row) {
  Bac_id <- row[["Bac"]]
  Euk_id <- row[["Euk"]]
  seas   <- row[["Season"]]
  ns     <- row[["NorthSouth"]]
  pairID <- row[["pair"]]

  # samples in this stratum
  sids <- meta2 %>%
    dplyr::filter(Season == seas, NorthSouth == ns) %>%
    dplyr::pull(SampleID)

  # intersect with available samples in matrices
  sids <- intersect(sids, colnames(X_bac))
  sids <- intersect(sids, colnames(X_euk))

  if (length(sids) == 0) return(NULL)

  # abundances
  bac_ab <- X_bac[Bac_id, sids, drop = TRUE]
  euk_ab <- X_euk[Euk_id, sids, drop = TRUE]

  pair_present <- as.integer(bac_ab > 0 & euk_ab > 0)

  tibble(
    pair      = pairID,
    Bac       = Bac_id,
    Euk       = Euk_id,
    SampleID  = sids,
    Season    = meta2$Season[match(sids, meta2$SampleID)],
    NorthSouth= meta2$NorthSouth[match(sids, meta2$SampleID)],
    Pair_present = pair_present
  )
}

# build the full sample-level table
sample_level <- purrr::map_dfr(
  seq_len(nrow(edges_pairs)),
  ~ get_pair_sample_presence(edges_pairs[.x, ])
)
message("matrices and metadata are ready!")
## Summaries per pair for Season and Latitude

pair_season <- sample_level %>%
  dplyr::group_by(pair, Season) %>%
  dplyr::summarise(
    n_present = sum(Pair_present),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = Season,
    values_from = n_present,
    values_fill = 0
  )

# Make sure we have columns named "spring" and "winter"
if (!"spring" %in% colnames(pair_season)) pair_season$spring <- 0L
if (!"winter" %in% colnames(pair_season)) pair_season$winter <- 0L


pair_lat <- sample_level %>%
  dplyr::group_by(pair, NorthSouth) %>%
  dplyr::summarise(
    n_present = sum(Pair_present),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = NorthSouth,
    values_from = n_present,
    values_fill = 0
  )

# ensure columns "North" and "South" exist
if (!"North" %in% colnames(pair_lat)) pair_lat$North <- 0L
if (!"South" %in% colnames(pair_lat)) pair_lat$South <- 0L

n_spring <- sum(meta2$Season == "spring")
n_winter <- sum(meta2$Season == "winter")

n_North  <- sum(meta2$NorthSouth == "North")
n_South  <- sum(meta2$NorthSouth == "South")

# Season enrichment per pair
pair_season <- pair_season %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    total_pair_season = spring + winter,
    p_season = if (total_pair_season > 0) {
      mat <- matrix(
        c(spring,
          winter,
          n_spring - spring,
          n_winter - winter),
        nrow = 2, byrow = TRUE
      )
      fisher.test(mat)$p.value
    } else {
      NA_real_
    },
    # log2 fold-change with +1 to avoid log2(0)
    log2FC_season = log2((spring + 1) / (winter + 1))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    q_season = p.adjust(p_season, method = "fdr")
  )

# Latitude enrichment per pair
pair_lat <- pair_lat %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    North = North,
    South = South,
    total_pair_lat = North + South,
    p_lat = if (total_pair_lat > 0) {
      mat <- matrix(
        c(North,
          South,
          n_North - North,
          n_South - South),
        nrow = 2, byrow = TRUE
      )
      fisher.test(mat)$p.value
    } else {
      NA_real_
    },
    log2FC_lat = log2((North + 1) / (South + 1))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    q_lat = p.adjust(p_lat, method = "fdr")
  )
message("Enrichment data ready")

## Combined enrichment table
enrich_all <- dplyr::full_join(
  pair_season %>% dplyr::select(pair, spring, winter, p_season, q_season, log2FC_season),
  pair_lat    %>% dplyr::select(pair, North,  South,  p_lat,    q_lat,    log2FC_lat),
  by = "pair"
)

# add Bac / Euk IDs and optional taxonomy
pair_ids <- edges_pairs %>%
  dplyr::distinct(pair, Bac, Euk)

enrich_all <- enrich_all %>%
  dplyr::left_join(pair_ids, by = "pair")

message("Enrichment data combined")

## Heatmap of enrichment
library(ggplot2)
library(forcats)

# select top N pairs by season or latitude significance
topN <- 250

top_pairs <- enrich_all %>%
  dplyr::arrange(q_season) %>%
  dplyr::slice_head(n = topN)

heat_df <- top_pairs %>%
  dplyr::select(pair, log2FC_season, log2FC_lat) %>%
  tidyr::pivot_longer(
    cols = c(log2FC_season, log2FC_lat),
    names_to = "contrast",
    values_to = "log2FC"
  )

heat_df$pair <- forcats::fct_reorder(heat_df$pair, heat_df$log2FC)

gg_heat <- ggplot(heat_df, aes(x = contrast, y = pair, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradient2(
    midpoint = 0,
    low  = "blue",
    mid  = "white",
    high = "red"
  ) +
  theme_bw(base_size = 10) +
  labs(x = "", y = "Pair (Bac||Euk)", fill = "log2FC") +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(OUTDIR, "enrichment_heatmap_top_pairs.pdf"),
       gg_heat, width = 6, height = 10)

message("Enrichment heatmap completed")


enrich_all_tax <- enrich_all %>%
  left_join(tax_bac, by = "Bac") %>%
  left_join(tax_euk, by = "Euk") %>%
  rename_with(~ paste0("Bac_", .), ends_with(".x")) %>%
  rename_with(~ paste0("Euk_", .), ends_with(".y"))

## Volcano plots => replace q per p (no correction)
top <- enrich_all_tax %>% arrange(p_season) %>% head(20)
gg_vol_season <- enrich_all %>%
  dplyr::filter(!is.na(p_season)) %>%
  ggplot(aes(x = log2FC_season, y = -log10(p_season))) +
  geom_point(alpha = 0.5) +
  aes(color = p_season < 0.05) +
  geom_text_repel(data = top, aes(label = paste(Bac_Taxonomy5.x, Euk_Taxonomy7.y)), max.overlaps = 40) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "Winter <- log2FC (Spring/Winter) -> Spring",
       y = "-log10(p_season)",
       title = "Seasonal enrichment")

top <- enrich_all_tax %>% arrange(p_lat) %>% head(20)
gg_vol_lat <- enrich_all %>%
  dplyr::filter(!is.na(p_lat)) %>%
  ggplot(aes(x = log2FC_lat, y = -log10(p_lat))) +
  geom_point(alpha = 0.5) +
  aes(color = p_lat < 0.05) +
  geom_text_repel(data = top, aes(label = paste(Bac_Taxonomy5.x, Euk_Taxonomy7.y)), max.overlaps = 40) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "South <- log2FC (North/South) -> North",
       y = "-log10(p_lat)",
       title = "Latitudinal enrichment")

volcano_combined <- gg_vol_season / gg_vol_lat  # patchwork stacking

ggsave(
  filename = file.path(OUTDIR, "volcano_plots_combined.pdf"),
  plot = gg_vol_season / gg_vol_lat,
  width = 8,
  height = 12
)

## collapsed Volcano plots
season_counts <- enrich_all_tax %>%
  dplyr::filter(!is.na(p_season)) %>%
  dplyr::group_by(spring, winter) %>%
  dplyr::summarise(
    n_pairs       = n(),
    n_pairs_log       = log10(n()),
    category = dplyr::if_else(
      n_pairs < 10,
      "Few (<10)",
      "Many (≥10)"),
    p_season      = first(p_season),
    log2FC_season = first(log2FC_season),
    .groups = "drop"
  ) %>%
  dplyr::mutate(mlog10p = -log10(p_season))

lat_counts <- enrich_all_tax %>%
  dplyr::filter(!is.na(p_lat)) %>%
  dplyr::group_by(North, South) %>%
  dplyr::summarise(
    n_pairs       = n(),
    n_pairs_log       = log10(n()),
    category = dplyr::if_else(
      n_pairs < 10,
      "Few (<10)",
      "Many (≥10)"
    ),
    p_lat      = first(p_lat),
    log2FC_lat = first(log2FC_lat),
    .groups = "drop"
  ) %>%
  dplyr::mutate(mlog10p = -log10(p_lat))

gg_vol_seasoncoll <- ggplot(season_counts, aes(x = log2FC_season, y = mlog10p, size = n_pairs_log)) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  ggplot2::scale_size_continuous(name = "Number of pairs", breaks = base::log10(c(1, 10, 100, 1000, 10000,100000)),   # log10 values
    labels = c("1", "10", "100", "1000", "10000",  "100000")) +
  ggplot2::geom_point(
    data = season_counts %>% dplyr::filter(n_pairs < 10),
    ggplot2::aes(),
    color = "grey60",
    alpha = 0.8
  ) +
  
  # Now draw common patterns with a gradient
  ggplot2::geom_point(
    data = season_counts %>% dplyr::filter(n_pairs >= 10),
    ggplot2::aes(color = n_pairs),
    alpha = 0.9
  ) +
  ggplot2::scale_color_gradient2(
    name = "Number of pairs",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 15000
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = "Winter <- log2FC (Spring/Winter) -> Spring",
    y = "-log10(p_season)",
    title = "Seasonal enrichment (collapsed volcano)"
  )
gg_vol_latcoll <- ggplot(lat_counts, aes(x = log2FC_lat, y = mlog10p, size = n_pairs_log)) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  ggplot2::scale_size_continuous(name = "Number of pairs", breaks = base::log10(c(1, 10, 100, 1000, 10000,100000)),   # log10 values
    labels = c("1", "10", "100", "1000", "10000",  "100000")) +
  ggplot2::geom_point(
    data = lat_counts %>% dplyr::filter(n_pairs < 10),
    ggplot2::aes(),
    color = "grey60",
    alpha = 0.8
  ) +
  
  # Now draw common patterns with a gradient
  ggplot2::geom_point(
    data = lat_counts %>% dplyr::filter(n_pairs >= 10),
    ggplot2::aes(color = n_pairs),
    alpha = 0.9
  ) +
  ggplot2::scale_color_gradient2(
    name = "Number of pairs",
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 15000
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = "South <- log2FC (North/South) -> North",
    y = "-log10(p_lat)",
    title = "Latitude enrichment (collapsed volcano)"
  )


ggsave(
  filename = file.path(OUTDIR, "volcano_coll_plots_combined.pdf"),
  plot = gg_vol_seasoncoll / gg_vol_latcoll,
  width = 8,
  height = 12
)

message("Volcano plots completed")


## Extract top season / latitude specific pairs
# Season-specific (FDR < 0.05)
sig_season <- enrich_all %>%
  dplyr::filter(q_season < 0.05) %>%
  dplyr::arrange(q_season)

# Latitude-specific
sig_lat <- enrich_all %>%
  dplyr::filter(q_lat < 0.05) %>%
  dplyr::arrange(q_lat)

# e.g., write them out:
readr::write_csv(sig_season, file.path(OUTDIR, "significant_pairs_season.csv"))
readr::write_csv(sig_lat,    file.path(OUTDIR, "significant_pairs_latitude.csv"))

message("Enrichment pipeline completed")

