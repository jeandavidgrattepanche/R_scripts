## with grazing data
setwd("/Volumes/JayDiii2SSD/Bioinfo_2025/R_2025")
# Output dir
OUTDIR    <- "plots/Nov25"


## MODEL A — Logistic enrichment for GrazingHL (High vs Low)
run_logistic_grazingHL <- function(df) {
  
  # Need both High and Low to fit model
  if (length(unique(df$GrazingHL)) < 2) {
    return(tibble::tibble(p = NA_real_, log2FC = NA_real_))
  }
  
  fit <- stats::glm(
    Pair_present ~ GrazingHL,
    data = df,
    family = stats::binomial()
  )
  
  coefs <- summary(fit)$coefficients
  
  # coefficient for High relative to Low
  p_val    <- coefs["GrazingHLHigh", "Pr(>|z|)"]
  log_odds <- coefs["GrazingHLHigh", "Estimate"]
  log2FC   <- log_odds / log(2)
  
  tibble::tibble(p = p_val, log2FC = log2FC)
}



## MODEL B — Logistic enrichment for Continuous Grazing
run_logistic_grazing_continuous <- function(df) {
  
  # Need variation
  if (stats::var(df$Grazing, na.rm = TRUE) == 0) {
    return(tibble::tibble(p = NA_real_, log2FC = NA_real_))
  }
  
  fit <- stats::glm(
    Pair_present ~ Grazing,
    data = df,
    family = stats::binomial()
  )
  
  coefs   <- summary(fit)$coefficients
  p_val   <- coefs["Grazing", "Pr(>|z|)"]
  log_odds <- coefs["Grazing", "Estimate"]
  log2FC   <- log_odds / log(2)   # log2 fold change per unit increase
  
  tibble::tibble(p = p_val, log2FC = log2FC)
}

#-------------- above are function -----------

## Add grazing info to sample_level
sample_level_g <- sample_level %>%
  dplyr::left_join(
    meta2 %>% dplyr::select(SampleID, Grazing, GrazingHL),
    by = "SampleID") %>%
  dplyr::mutate(
      GrazingHL = factor(GrazingHL, levels = c("Low", "High"))  # ensure Low is reference
  )
## compute grazing HL
grazingHL_results <- sample_level_g %>%
  dplyr::group_by(pair) %>%
  dplyr::group_modify(~ run_logistic_grazingHL(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(q = stats::p.adjust(p, method = "fdr"))

gg_vol_grazHL <- ggplot2::ggplot(
  grazingHL_results,
  ggplot2::aes(x = log2FC, y = -base::log10(p))
) +
  ggplot2::geom_point(alpha = 0.5) +
  ggplot2::geom_hline(yintercept = -base::log10(0.05), linetype = "dashed", color = "red") +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = "log2FC (High vs Low grazing)",
    y = "-log10(p)",
    title = "GrazingHL enrichment (High/Low)"
  )


## compute grazing continuous
grazing_cont_results <- sample_level_g %>%
  dplyr::group_by(pair) %>%
  dplyr::group_modify(~ run_logistic_grazing_continuous(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(q = stats::p.adjust(p, method = "fdr"))

gg_vol_grazcont <-ggplot2::ggplot(
  grazing_cont_results,
  ggplot2::aes(x = log2FC, y = -base::log10(p))
) +
  ggplot2::geom_point(alpha = 0.5) +
  ggplot2::geom_hline(yintercept = -base::log10(0.05), linetype = "dashed", color = "red") +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = "log2FC (per unit increase in grazing)",
    y = "-log10(p)",
    title = "Continuous grazing enrichment"
  )

ggsave(
  filename = file.path(OUTDIR, "volcano_grazing_plots_combined.pdf"),
  plot = gg_vol_grazHL / gg_vol_grazcont,
  width = 8,
  height = 12
)

message("Volcano plots completed")


