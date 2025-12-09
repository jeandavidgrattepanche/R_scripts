# reproducible R/phyloseq workflow to discover Bacteria↔Eukaryote pairs and test your three hypotheses:
# seasonality (winter vs spring)
# latitude (North vs South Antarctic Peninsula)
# effect of nanoflagellate grazing (rate/covariate)
# 
# R pipeline that assumes you already have two phyloseq objects (or tables you can import into them):
#     ps_bac — bacterial OTUs/ASVs + taxonomy + sample_data
#     ps_euk — eukaryote OTUs/ASVs + taxonomy + sample_data
# sample metadata in both has: 
#     Season ∈ {Winter, Spring}, 
#     Region ∈ {North, South}, 
#     Grazing (numeric; “bacteria consumed by nanoflagellates” per sample)
# You’ll get:
#     signed cross-domain correlations per stratum (season×region and grazing splits)
#     optional partial correlations controlling for Grazing
#     FDR control, edge lists (CSV)
#     heatmaps of taxon–taxon interaction intensity
#     quick stats for seasonal/latitudinal enrichment 


setwd("/Volumes/JayDiii2SSD/Bioinfo_2025/R_2025")
# Output dir
OUTDIR    <- "plots/Nov25"
dir.create(OUTDIR, showWarnings = FALSE)
SEED         <- 1

pkgs <- c(
  "phyloseq","dplyr","tidyr","tibble","purrr","stringr","phyloseqCompanion","plyr","ape","phytools",
  "microbiome","Matrix","psych","ppcor","ggplot2","ggrepel","pheatmap","forcats","future.apply","readr","patchwork"
)
invisible(lapply(pkgs, function(p){ if (!requireNamespace(p, quietly=TRUE)) { install.packages(p)}}))
invisible(lapply(pkgs, function(p) {suppressPackageStartupMessages(suppressMessages(library(p, character.only=TRUE)))}))

set.seed(SEED)


#=>  run Load_data_nov25.R first
message("Load data 303 ")
source('Load_data_nov25.r', chdir = TRUE)

# can module size here (maybe no micro?) ps_micro_r
 ps_euk <- merge_phyloseq(ps_pico_r, ps_nano_r)


# keep only shared samples
common_samples <- intersect(sample_names(ps_bac_r), sample_names(ps_euk))
ps_bac <- prune_samples(common_samples, ps_bac_r)
ps_euk <- prune_samples(common_samples, ps_euk)

# ensure key metadata exist => water_temp will be replace by Grazing when script is working properly
stopifnot(all(c("Season","NorthSouth","latitude") %in% colnames(sample_data(ps_bac))))
stopifnot(all(c("Season","NorthSouth","latitude") %in% colnames(sample_data(ps_euk))))

message("Data loaded and verify for Season, latitude and true latitude (to be replace by grazing later)")

# (optional) taxonomic agglomeration to desired ranks => off to keep species level
# change "Order" / "Class" to what makes sense for your data
# ps_bac_ord <- tax_glom(ps_bac, taxrank = "Taxonomy6", NArm = TRUE)
# ps_euk_cls <- tax_glom(ps_euk, taxrank = "Taxonomy6", NArm = TRUE)

# filter rare features
# ps_bac_ordr <- prune_taxa(taxa_sums(ps_bac_ord) > 10, ps_bac_ord)
# ps_euk_clsr <- prune_taxa(taxa_sums(ps_euk_cls) > 10, ps_euk_cls)
ps_bac_r <- prune_taxa(taxa_sums(ps_bac) > 10, ps_bac)
ps_euk_r <- prune_taxa(taxa_sums(ps_euk) > 10, ps_euk)

#filter by taxa
keep_euk_tax4 <- c(
  "p:Stramenopiles-Gyrista",
  "p:Stramenopiles-Bigyra",
  "p:Stramenopiles-Stramenopiles_X",
  "p:Stramenopiles:plas-Gyrista:plas",
  "p:Alveolata-Ciliophora",
  "p:Alveolata-Dinoflagellata",
  "p:Alveolata-Alveolata_X",
  "p:Cryptophyta-Cryptophyta_X",
  "p:Cryptophyta:nucl-Cryptophyta_X:nucl",
  "p:Haptophyta-Haptophyta_X",
  "p:Picozoa-Picozoa_X",
  "p:Telonemia-Telonemia_X",
  "p:Rhizaria-Cercozoa",
  "p:Rhizaria-Radiolaria",
  "p:Opisthokonta-Choanoflagellata",
  "p:Discosea-Discosea_X",
  "p:Tubulinea-Tubulinea_X",
  "p:Centroplasthelida-Centroplasthelida_X",
  "p:Prasinodermophyta-Prasinodermophyta_X"
)
keep_tax5 <- c(
  # Diatoms
  "c:Mediophyceae", "c:Coscinodiscophyceae", "c:Bacillariophyceae",
  "c:Mediophyceae:plas", "c:Diatomeae_X",

  # Chrysophytes, Dictyochophytes, Pelagophytes, MOCH
  "c:Chrysophyceae", "c:Dictyochophyceae", "c:Pelagophyceae",
  "c:Bolidophyceae", "c:Opalozoa", "c:MOCH-1", "c:MOCH-2", "c:Sagenista",

  # Cryptophytes
  "c:Cryptophyceae", "c:Cryptophyceae:nucl",

  # Haptophytes
  "c:Prymnesiophyceae", "c:Haptophyta_XX",

  # Dinoflagellates
  "c:Dinophyceae", "c:Dinophyta_X", "c:Syndiniales",

  # Green picoeuks
  "c:Mamiellophyceae", "c:Trebouxiophyceae",
  "c:Prasinodermophyceae", "c:Pyramimonadophyceae",

  # Picozoa
  "c:Picozoa_XX", "c:Picomonadea",

  # Telonemia
  "c:Telonemia_XX",

  # Choanoflagellates
  "c:Choanoflagellatea",

  # Rhizaria (predators)
  "c:Acantharea", "c:Polycystinea", "c:Cercozoa_X",
  "c:Filosa-Thecofilosea", "c:Filosa-Imbricatea",
  "c:Filosa-Sarcomonadea", "c:Filosa-Granofilosea",
  "c:Phaeodarea", "c:Endomyxa-Ascetosporea",
  "c:Rotosphaerida_X", "c:RAD-B", "c:RAD-C",

  # Amoebozoa
  "c:Discosea_XX", "c:Centramoebia",

  # environmental / unclassified Stramenopiles
  "c:Gyrista_X", "c:Stramenopiles_XX",
  "c:CONThreeP", "c:CONTH_7",
  "c:Stramenopiles_X-Group-2", "c:Stramenopiles_X-Group-9"
)

keep_tax5_no_diatoms <- c(
  # Chrysophytes, Dictyochophytes, Pelagophytes, MOCH
  "c:Chrysophyceae", "c:Dictyochophyceae", "c:Pelagophyceae",
  "c:Bolidophyceae", "c:Opalozoa", "c:MOCH-1", "c:MOCH-2", "c:Sagenista",

  # Cryptophytes
  "c:Cryptophyceae", "c:Cryptophyceae:nucl",

  # Haptophytes
  "c:Prymnesiophyceae", "c:Haptophyta_XX",

  # Dinoflagellates
  "c:Dinophyceae", "c:Dinophyta_X", "c:Syndiniales",

  # Green picoeuks
  "c:Mamiellophyceae", "c:Trebouxiophyceae",
  "c:Prasinodermophyceae", "c:Pyramimonadophyceae",

  # Picozoa
  "c:Picozoa_XX", "c:Picomonadea",

  # Telonemia
  "c:Telonemia_XX",

  # Choanoflagellates
  "c:Choanoflagellatea",

  # Rhizaria (predators)
  "c:Acantharea", "c:Polycystinea", "c:Cercozoa_X",
  "c:Filosa-Thecofilosea", "c:Filosa-Imbricatea",
  "c:Filosa-Sarcomonadea", "c:Filosa-Granofilosea",
  "c:Phaeodarea", "c:Endomyxa-Ascetosporea",
  "c:Rotosphaerida_X", "c:RAD-B", "c:RAD-C",

  # Amoebozoa
  "c:Discosea_XX", "c:Centramoebia",

  # environmental / unclassified Stramenopiles
  "c:Gyrista_X", "c:Stramenopiles_XX",
  "c:CONThreeP", "c:CONTH_7",
  "c:Stramenopiles_X-Group-2", "c:Stramenopiles_X-Group-9"
)
keep_tax6 <- c(
  # Heterotrophic micrograzers & mixotrophs
  "o:Oligotrichida", "o:Choreotrichida", "o:Choreotrichida-Tintinnina",
  "o:Oligohymenophorea_X", "o:Hypotrichia", "o:Scuticociliatia_1",
  "o:Haptoria_4", "o:Haptoria_6", "o:Haptoria_7",
  "o:Euplotia", "o:Prostomatea", "o:Nassophorea_X", "o:Cyrtophoria_1",
  "o:Cyrtophoria_2", "o:Cyrtophoria_8",

  # Dinoflagellates (mixotrophs, predators) "o:Gymnodiniales" <- removed for now
  "o:Peridiniales", "o:Torodiniales", "o:Prorocentrales",
  "o:Dinophysiales", "o:Gonyaulacales", "o:Suessiales",
  "o:Dino-Group-I", "o:Dino-Group-II", "o:Dino-Group-III",
  "o:Dino-Group-IV", "o:Dino-Group-V", "o:Dinophyta_XX",

  # Cryptophytes
  "o:Cryptomonadales", "o:Cryptomonadales:nucl", "o:Pyrenomonadales",

  # Chrysophytes / Dictyochophytes / Pelagophytes
  "o:Chrysophyceae_Clade-EC2H", "o:Chrysophyceae_Clade-EC2I",
  "o:Chrysophyceae_X", "o:Sarcinochrysidales",
  "o:Dictyochophyceae_X", "o:Sagenista_X",
  "o:Pelagomonadales", "o:Ochromonadales", "o:Paraphysomonadales",
  "o:Dolichomastigales",

  # Haptophytes (often mixotrophic)
  "o:Prymnesiales", "o:Phaeocystales",
  "o:Prymnesiophyceae_Clade_E", "o:Prymnesiophyceae_Clade_D",
  "o:Prymnesiophyceae_Clade_F",
  "o:Calcihaptophycidae",

  # Picozoa
  "o:Picozoa_XXX", "o:Picomonadida",

  # Prasinophytes / pico green algae
  "o:Mamiellales", "o:Pyramimonadales", "o:Prasinodermales",
  "o:Chlamydomonadales",

  # Rhizaria (major heterotrophic grazers)
  "o:Acantharea_4", "o:Acantharea_3", "o:Acantharea_B", "o:Acantharea_C", "o:Acantharea_F",
  "o:Polycystinea", "o:Spumellaria", "o:Nassellaria",
  "o:Rhabdonematales", "o:Phaeogromida",
  "o:Cercozoa_XX", "o:Thaumatomonadida", "o:Euglyphida",
  "o:Cercomonadida", "o:Glissomonadida",
  "o:Filosa-Imbricatea_X", "o:Filosa-Thecofilosea_X",
  "o:Filosa-Granofilosea_X",
  "o:Endomyxa-Ascetosporea", "o:Pseudoperkinsidae",

  # Telonemia
  "o:Telonemia_XXX",

  # Labyrinthulomycetes (osmotrophs, but ecologically linked)
  "o:Labyrinthulomycetes", "o:Developea_X",

  # Opalozoa
  "o:Opalozoa_X",

  # Choanoflagellates (key grazers)
  "o:Craspedida", "o:Acanthoecida", "o:Choanoflagellatea_X",

  # Amoebozoa
  "o:Vannellida", "o:Centramoebia", "o:Leptomyxida",

  # Environmental Stramenopiles (good to keep)
  "o:Gyrista_XX", "o:Stramenopiles_XXX",
  "o:Stramenopiles_X-Group-2_X", "o:Stramenopiles_X-Group-9_X",
  "o:CONThreeP_X", "o:CONTH_7_X", "o:OLIGO5", "o:Pseudophyllomitidae",

  # Parmales (tiny algae; keep if studying mixotrophy)
  "o:Parmales",

  # Picoeuk taxa with ecological relevance
  "o:Nanomonadea", "o:Marimonadida", "o:Picomonadida"
)

bad_tax4_exact <- c(
  "Chloroplast",
  "Chlorophyta",
  "Dinoflagellata",
  "Ochrophyta",
  "Phaeocystis",
  "Ciliophora",
  "Cercozoa",
  "Retaria",
  "Opisthokonta",
  "Telonema_antarcticum",
  "Prymnesiales",
  "Coccolithales",
  "Cryptophyta_sp._SL64/78sp",
  "Geminigera",
  "MAST-1",
  "MAST-3",
  "MAST-7"
)

bad_tax4_patterns <- c(
  "Chloroplast",
  "Chlorophyta",
  "Dinoflagellata",
  "Ochrophyta",
  "Phaeocystis",
  "Ciliophora",
  "Cercozoa",
  "Retaria",
  "Opisthokonta",
  "Telonema",
  "Prymnesiales",
  "Coccolithales",
  "Cryptophyta",
  "MAST-",
  "Geminigera"
)


# modify taxa of interest => but want the species level name for both bact and euk
ps_euk_filt <- subset_taxa(ps_euk_r, Taxonomy6 %in% keep_tax6)
# ps_bac_filt <- subset_taxa(ps_bac_ordr, !(Taxonomy4 %in% bad_tax4_exact))
ps_bac_filt <- subset_taxa(
  ps_bac_r,
  {
    tx4 <- as.character(Taxonomy4)
    # TRUE if Taxonomy4 contains any "bad" pattern
    is_bad <- Reduce(`|`, lapply(bad_tax4_patterns, function(p) grepl(p, tx4)))
    !is_bad  # keep only non-bad
  }
)
message("Data filtered Bac => tax4 and Euk => tax6")

# 
# ntaxa(ps_euk_cls)
# ntaxa(ps_euk_filt)
# ntaxa(ps_bac_ordr)
# ntaxa(ps_bac_filt)

### CLR transform (composition-aware)
clr_transform <- function(ps, pseudocount = 1) {
  mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) mat <- mat else mat <- t(mat)
  mat <- mat + pseudocount
  mat <- microbiome::transform(mat, "clr")
  # return features in rows, samples in cols
  mat
}

X_bac <- clr_transform(ps_bac_filt)  # features x samples
X_euk <- clr_transform(ps_euk_filt)  # features x samples
message("CLR transformation done")


meta <- as(sample_data(ps_bac_filt), "data.frame") %>% 
  rownames_to_column("SampleID")

## helper: cross-domain correlation in a stratum mode fast

plan(multisession, workers = parallel::detectCores() - 1)

corr_bac_euk_fast <- function(X_bac, X_euk, sample_ids,
                              method    = "spearman",
                              fdr_alpha = 0.05,
                              min_abs_r = 0.3) {
  # intersect sample ids with matrices
  sid <- intersect(sample_ids, intersect(colnames(X_bac), colnames(X_euk)))
  if (length(sid) < 3) {
    # not enough samples: return empty *tibble*
    return(tibble(
      Bac = character(0),
      Euk = character(0),
      rho = numeric(0),
      p   = numeric(0),
      q   = numeric(0)
    ))
  }

  B <- X_bac[, sid, drop = FALSE]
  E <- X_euk[, sid, drop = FALSE]
  n <- length(sid)

  res_list <- future_lapply(seq_len(nrow(B)), function(i) {
    r_vec <- apply(E, 1, function(e) suppressWarnings(cor(B[i, ], e, method = method)))
    # approximate p from r
    # guard against |r| = 1
    r_vec[abs(r_vec) >= 1] <- sign(r_vec[abs(r_vec) >= 1]) * (1 - 1e-10)
    t_vec <- r_vec * sqrt(n - 2) / sqrt(1 - r_vec^2)
    p_vec <- 2 * pt(-abs(t_vec), df = n - 2)

    tibble(
      Bac = rownames(B)[i],
      Euk = rownames(E),
      rho = as.numeric(r_vec),
      p   = as.numeric(p_vec)
    )
  })

  df <- bind_rows(res_list) %>%
    dplyr::mutate(q = p.adjust(p, "BH")) %>%
    filter(!is.na(rho), !is.na(q),
           abs(rho) >= min_abs_r,
           q <= fdr_alpha)

  df
}

## run per stratum (Season × Region) => modify water_temp by actual grazing data
meta2 <- meta %>%
  dplyr::mutate(Season = factor(Season),
         NorthSouth = factor(NorthSouth),
         Grazing = as.numeric(latitude),
         GrazingHL = if_else(Grazing >= median(Grazing, na.rm=TRUE), "High","Low"))

# Season × Region
strata_SR <- meta2 %>%
  dplyr::group_by(Season, NorthSouth) %>%
  dplyr::summarize(n = dplyr::n(), Samples = list(SampleID), .groups="drop")

edges_SR <- strata_SR %>%
  dplyr::mutate(
    edges = purrr::map(
      Samples,
      ~ corr_bac_euk_fast(
          X_bac, X_euk,
          sample_ids = .x,
          min_abs_r  = 0.2,  # same as your test
          fdr_alpha  = 0.1
        )
    )
  ) %>%
  dplyr::select(Season, NorthSouth, n, edges) %>%
  tidyr::unnest(cols = edges)

message("SR network data ready")

# Grazing split
strata_G <- meta2 %>%
  dplyr::group_by(GrazingHL) %>%
  dplyr::summarize(n = dplyr::n(), Samples = list(SampleID), .groups="drop")

edges_G <- strata_G %>%
#   dplyr::mutate(res = purrr::map(Samples, ~ corr_bac_euk_fast(X_bac, X_euk, .x))) %>%
#   dplyr::mutate(edges = purrr::map(res, "edges")) %>%
  dplyr::mutate(edges = purrr::map(Samples, ~ corr_bac_euk_fast(X_bac, X_euk, .x, min_abs_r  = 0.2, fdr_alpha  = 0.1))) %>%
  dplyr::select(GrazingHL, n, edges) %>%
  tidyr::unnest(cols = edges)

message("G network data ready")


## annotate edges with taxonomy and summarize by rank
tax_bac <- as(tax_table(ps_bac_filt), "matrix") %>% as.data.frame() %>% rownames_to_column("Bac")
tax_euk <- as(tax_table(ps_euk_filt), "matrix") %>% as.data.frame() %>% rownames_to_column("Euk")

annotate_edges <- function(edges_df) {
  edges_df %>%
    left_join(tax_bac, by="Bac") %>%
    left_join(tax_euk, by="Euk", suffix=c("_Bac","_Euk"))
}

edges_SR_annot <- annotate_edges(edges_SR)
edges_G_annot  <- annotate_edges(edges_G)

message("edges annotation done")


# heatmap prep: aggregate by chosen ranks (here Order_Bac × Class_Euk)
agg_SR <- edges_SR_annot %>%
  dplyr::group_by(Season, NorthSouth, species_Bac, species_Euk) %>%
  dplyr::summarize(n_edges = dplyr::n(),
            rho_mean = mean(rho, na.rm=TRUE),
            .groups="drop")

## heatmaps (per stratum)
plot_heatmap_sr <- function(df, fill = "n_edges") {
  # one heatmap per Season×Region
  splits <- df %>% dplyr::group_split(Season, NorthSouth)
  plots <- lapply(splits, function(d) {
    this_season <- unique(d$Season)
    this_region <- unique(d$NorthSouth)
    message("Plotting: ", this_season, " - ", this_region)
    M <- d %>% 
        dplyr::select(species_Bac, species_Euk, value = !!rlang::sym(fill)) %>%
        tidyr::pivot_wider(names_from = species_Euk, values_from = value, values_fill = list(value =0)) %>%
        tibble::column_to_rownames("species_Bac") %>% as.matrix()
    pheatmap::pheatmap(M, cluster_rows=TRUE, cluster_cols=TRUE, main = paste0(unique(d$Season), " - ", unique(d$NorthSouth)), fontsize = 9, border_color = "grey90")
  })
  invisible(plots)
}

# example: number of significant edges
pdf(file.path(OUTDIR,"heatmap_network.pdf"))
plot_heatmap_sr(agg_SR, fill = "n_edges")
dev.off()



## remove meaningless info here below 20%
cutoff_percent <- 0.01

# compute cutoff per S and R
agg_SR_cut <- agg_SR %>%
  dplyr::group_by(Season, NorthSouth) %>%
  dplyr::mutate(
        cutoff_value= quantile(n_edges, probs = 1 - cutoff_percent, na.rm= TRUE),
        keep = n_edges >= cutoff_value
  ) %>%
  dplyr::ungroup()%>%
  dplyr::filter(keep)

plot_heatmap_sr_percent <- function(df, fill = "n_edges") {
  # one heatmap per Season×Region
  splits <- df %>% dplyr::group_split(Season, NorthSouth)
  plots <- lapply(splits, function(d) {
    this_season <- unique(d$Season)
    this_region <- unique(d$NorthSouth)
    message("Plotting: ", this_season, " - ", this_region)
    mat <- d %>% 
        dplyr::select(species_Bac, species_Euk, value = !!rlang::sym(fill)) %>%
        tidyr::pivot_wider(names_from = species_Euk, values_from = value, values_fill = list(value =0)) %>%
        tibble::column_to_rownames("species_Bac") %>% as.matrix()
    pheatmap::pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, main = paste0(unique(d$Season), " - ", unique(d$NorthSouth)), fontsize = 9, border_color = "grey90", color = viridisLite::viridis(50))
  })
  invisible(plots)
}

# example: number of significant edges
pdf(file.path(OUTDIR,"heatmap_network_percent.pdf"))
plot_heatmap_sr_percent(agg_SR_cut, fill = "n_edges")
dev.off()


## keep only top 20 interaction
topN <- 30

# compute topN per S and R
agg_SR_topN <- agg_SR %>%
  dplyr::group_by(Season, NorthSouth) %>%
  dplyr::arrange(dplyr::desc(n_edges)) %>%
  dplyr::slice_head(n = topN)%>%
  dplyr::ungroup()
  
plot_heatmap_sr_topN <- function(df, fill = "n_edges") {
  # one heatmap per Season×Region
  splits <- df %>% dplyr::group_split(Season, NorthSouth)
  plots <- lapply(splits, function(d) {
    this_season <- unique(d$Season)
    this_region <- unique(d$NorthSouth)
    message("Plotting: ", this_season, " - ", this_region)
    mat <- d %>% 
        dplyr::select(species_Bac, species_Euk, value = !!rlang::sym(fill)) %>%
        tidyr::pivot_wider(names_from = species_Euk, values_from = value, values_fill = list(value =0)) %>%
        tibble::column_to_rownames("species_Bac") %>% as.matrix()
    pheatmap::pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE, main = paste0(unique(d$Season), " - ", unique(d$NorthSouth)), fontsize = 9, border_color = "grey90", color = viridisLite::viridis(80))
  })
  invisible(plots)
}

# example: number of significant edges
pdf(file.path(OUTDIR,"heatmap_network_topN.pdf"))
plot_heatmap_sr_topN(agg_SR_topN, fill = "n_edges")
dev.off()

# # check for data for each pear
# split_SR <- edges_SR_annot %>% dplyr::group_split(Season, NorthSouth)
# 
# names(split_SR) <- edges_SR_annot %>% 
#   dplyr::distinct(Season, NorthSouth) %>% 
#   dplyr::mutate(name = paste(Season, NorthSouth, sep = "_")) %>% 
#   dplyr::pull(name)
# 
# lapply(split_SR, nrow)

## save edge lists
readr::write_csv(edges_SR_annot, file.path(OUTDIR,"edges_significant_by_SeasonRegion_relax.csv"))
readr::write_csv(edges_G_annot,  file.path(OUTDIR,"edges_significant_by_GrazingHL.csv"))

