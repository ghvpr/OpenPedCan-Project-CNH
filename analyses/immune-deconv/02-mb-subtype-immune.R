############################################
# 02-mb-subtype-immune.R
#
# Purpose:
# Compare immune cell distributions across
# medulloblastoma molecular subtypes
# (WNT, SHH, Group3, Group4) using
# immune deconvolution results and
# histology annotations.
############################################

# ---- Load libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

# ---- Define input paths ----
# Paths are relative to analyses/immune-deconv/
hist_path   <- "/home/rstudio/project/data/histologies.tsv"
deconv_path <- "/home/rstudio/project/analyses/immune-deconv/results/quantiseq_output.rds"

# ---- Read input files ----
hist   <- read_tsv(hist_path, show_col_types = FALSE)
deconv <- readRDS(deconv_path)

# ---- Sanity checks ----
cat("Histologies table rows:", nrow(hist), "\n")
cat("Deconvolution table rows:", nrow(deconv), "\n")

cat("Immune deconvolution methods detected:\n")
print(unique(deconv$method))

mb_rows <- sum(grepl("medullo", deconv$cancer_group, ignore.case = TRUE))
cat("Rows labeled as medulloblastoma in deconvolution output:", mb_rows, "\n")


# 1) Keep medulloblastoma samples using short_histology
# 2) Map molecular_subtype -> WNT/SHH/Group3/Group4
# 3) Keep only the join key + subtype label, one row per biospecimen

hist_mb_subtype <- hist %>%
  # medulloblastoma-only
  filter(str_detect(short_histology, regex("medullo", ignore_case = TRUE))) %>%
  # clean subtype labels
  mutate(mb_subtype = case_when(
    str_detect(molecular_subtype, regex("WNT",    ignore_case = TRUE)) ~ "WNT",
    str_detect(molecular_subtype, regex("SHH",    ignore_case = TRUE)) ~ "SHH",
    str_detect(molecular_subtype, regex("Group3", ignore_case = TRUE)) ~ "Group3",
    str_detect(molecular_subtype, regex("Group4", ignore_case = TRUE)) ~ "Group4",
    TRUE ~ NA_character_
  )) %>%
  # drop unclassified / other
  filter(!is.na(mb_subtype)) %>%
  # keep only what we need for the join
  select(Kids_First_Biospecimen_ID, mb_subtype) %>%
  distinct()


cat("MB subtype annotation rows:", nrow(hist_mb_subtype), "\n")
print(table(hist_mb_subtype$mb_subtype))
cat("Duplicate biospecimen IDs:", anyDuplicated(hist_mb_subtype$Kids_First_Biospecimen_ID), "\n")

deconv_filt <- deconv %>%
  filter(method == "quantiseq") %>%
  filter(sample_type == "Tumor")

cat("Deconv rows after filtering (quantiseq + Tumor):", nrow(deconv_filt), "\n")

mb_immune <- deconv_filt %>%
  inner_join(hist_mb_subtype, by = "Kids_First_Biospecimen_ID")

cat("Rows after join (MB only):", nrow(mb_immune), "\n")
print(table(mb_immune$mb_subtype))

# (3) Output directories
out_dir <- "results/mb-subtype-immune"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# (4) Summarize mean fraction per subtype & immune cell type
mb_summary <- mb_immune %>%
  group_by(mb_subtype, cell_type) %>%
  summarize(
    n = n(),
    mean_fraction = mean(fraction, na.rm = TRUE),
    median_fraction = median(fraction, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(mb_summary, file.path(out_dir, "mb_subtype_quantiseq_summary.tsv"))


p_stacked <- mb_summary %>%
  filter(cell_type !="uncharacterized cell") %>%
  group_by(mb_subtype) %>%
  mutate(prop_of_mean = mean_fraction / sum(mean_fraction)) %>%
  ungroup() %>%
  ggplot(aes(x = mb_subtype, y = prop_of_mean, fill = cell_type)) +
  geom_col() +
  labs(
    title = "Mean immune cell composition by medulloblastoma subtype (quanTIseq), excluding uncharacterized cells",
    x = "Medulloblastoma subtype",
    y = "Proportion of mean immune fraction",
    fill = "Immune cell type"
  ) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "mb_subtype_mean_composition_stacked_characterized.png"),
  plot = p_stacked,
  width = 10, height = 6, dpi = 300
)
