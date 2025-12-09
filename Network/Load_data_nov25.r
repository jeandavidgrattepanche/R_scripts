# ============================================================
# load_phyloseq_nov25.R
# Loads B303 bacterial + eukaryotic datasets into phyloseq objects
# Returns a list of phyloseq objects.
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
OUTDIR    <- "plots/Nov25"
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
pkgs <- c(
  "phyloseq","phyloseqCompanion","phytools","tidyverse"
)
invisible(lapply(pkgs, function(p){ if (!requireNamespace(p, quietly=TRUE)) { install.packages(p)}}))
invisible(lapply(pkgs, function(p) {suppressPackageStartupMessages(suppressMessages(library(p, character.only=TRUE)))}))

set.seed(SEED)

# ============================================================
# ---- 2) Helpers --------------------------------------------

.clean_ps <- function(ps) {
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
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
