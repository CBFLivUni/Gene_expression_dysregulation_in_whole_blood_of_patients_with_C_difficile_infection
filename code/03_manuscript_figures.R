## Set up ####
### Packages used ####
library(tidytable)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(forcats)
library(oligo)
library(conflicted)
library(gridExtra)
library(janitor)
library(fgsea)

conflicts_prefer(tidytable::rename,
                 tidytable::filter)

### Plot theme ####
# Colour palette
okabe_ito <- as.vector(palette.colors())

# Custom theme_bw to remove dark grey colours
theme_bw_black <- function(base_size = 14,
                       base_family = "",
                       base_line_size = base_size / 22,
                       base_rect_size = base_size / 22)
{
  theme_grey(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) +
    theme(
      line = element_line(color = "black"),
      text = element_text(color = "black"),
      #   axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "white",
                                      colour = NA),
      panel.border = element_rect(fill = NA,
                                  colour = "black"),
      panel.grid = element_line(colour = "grey92"),
      panel.grid.minor = element_line(linewidth = rel(0.5)),
      strip.background = element_rect(fill = "white",
                                      colour = "black"),
      complete = FALSE
    )
}

## Figure 2 ####
### Volcano plot functions ####
# Volcano plot function
volcano <- function(df, top20 = NULL) {
  plt <- df |> 
    mutate(highlight = case_when(
      (pval_under_05 & abs(logFC) > 0.5) ~ "C",
      pval_under_05 ~ "B",
      TRUE ~ "A"
    )) |> 
    ggplot(aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(colour = highlight),
               alpha = 0.8,
               # shape = 21,
               #   colour = "white",
               #   stroke = 0.01
    ) +
    scale_colour_manual(values = c(as.vector(palette.colors())[c(1,3)], "dodgerblue4")) +
    theme_bw_black() +
    theme(legend.position = "none") +
    labs(x = expression(Log[2]~fold~change),
         y = expression("−"*Log[10]*"(p-value)")
    )
  
  if (!is.null(top20)) {
    plt <- plt + geom_text_repel(data = top20 |> 
                                   filter(adj.P.Val < 0.05),
                                 aes(label = SYMBOL),
                                 size = 2,
                                 max.overlaps = 20)
  }
  return(plt)
}

# Function to calculate the number of upregulated and downregulated DE genes
# and format as a table
de_table <- function(df, heading) {
  # Create a new data frame with the desired columns
  de_count <- data.frame(
    direction = c("↑", "↓"),
    number = c(
      sum(df$adj.P.Val < 0.05 & df$logFC > 0),
      sum(df$adj.P.Val < 0.05 & df$logFC < 0)
    )
  )
  # Create one-column data frame with heading
  de_table <- data.frame(
    count_string = paste(de_count$direction, de_count$number)
  )
  
  colnames(de_table) <- paste(heading)
  
  return(de_table)
}

### CDI vs. HC ####
# Import data
tbl_cdi_hc_bsln <- fread("output/tables/de_stats_cdi_vs_hc.csv")

# Identify top20 upregulated  genes for point labelling
top20_cdi_hc_bsln <- tbl_cdi_hc_bsln |> 
  arrange(desc(logFC)) |> 
  slice_head(20)

# Plot
vlcno_cdi_hc <- volcano(tbl_cdi_hc_bsln, top20_cdi_hc_bsln) +
  annotation_custom(tableGrob(de_table(tbl_cdi_hc_bsln, "CDI vs. HC"),
                              theme = ttheme_minimal(
                                core = list(
                                  fg_params = list(hjust = 0.5))),
                              rows = NULL),
                    xmin = 2.7, xmax = 3.5,
                    ymin = 2.1, ymax = 6)

# Save plot
save_plot("output/figures/volcano_cdi_vs_hc.pdf",
          vlcno_cdi_hc,
          base_asp = 1.2, 
          ncol = 1,
          device = cairo_pdf)

### GDH+ vs. HC ####
# Import data
tbl_gdh_hc_bsln <- fread("output/tables/de_stats_gdh_vs_hc.csv")

# Identify top20 upregulated  genes for point labelling
top20_gdh_hc_bsln <- tbl_gdh_hc_bsln |> 
  arrange(desc(logFC)) |> 
  slice_head(20)

# Plot
vlcno_gdh_hc <- volcano(tbl_gdh_hc_bsln, top20_gdh_hc_bsln) +
  annotation_custom(tableGrob(de_table(tbl_gdh_hc_bsln, "GDH vs.\nHC"),
                              theme = ttheme_minimal(
                                core = list(
                                  fg_params = list(hjust = 0.5))),
                              rows = NULL),
                    xmin = 2.1, xmax = 3,
                    ymin = 0.6, ymax = 3.5)

# Save plot
save_plot("output/figures/volcano_gdh_vs_hc.pdf",
          vlcno_gdh_hc,
          base_asp = 1.2, 
          ncol = 1,
          device = cairo_pdf)

### IBD vs. HC ####
# Import data
tbl_ibd_hc_bsln <- fread("output/tables/de_stats_ibd_vs_hc.csv")

# Identify top20 upregulated  genes for point labelling
top20_ibd_hc_bsln <- tbl_ibd_hc_bsln |> 
  arrange(desc(logFC)) |> 
  slice_head(20)

# Plot
vlcno_ibd_hc <- volcano(tbl_ibd_hc_bsln, top20_ibd_hc_bsln) +
  annotation_custom(tableGrob(de_table(tbl_ibd_hc_bsln, "IBD vs. HC"),
                              theme = ttheme_minimal(
                                core = list(
                                  fg_params = list(hjust = 0.5))),
                              rows = NULL),
                    xmin = 3.2, xmax = 4.5,
                    ymin = 2.1, ymax = 6)

# Save plot
save_plot("output/figures/volcano_ibd_vs_hc.pdf",
          vlcno_ibd_hc,
          base_asp = 1.2, 
          ncol = 1,
          device = cairo_pdf)

### DC vs. HC ####
# Import data
tbl_dc_hc_bsln <- fread("output/tables/de_stats_dc_vs_hc.csv")

# Identify top20 upregulated genes for point labelling
top20_dc_hc_bsln <- tbl_dc_hc_bsln |> 
  arrange(desc(logFC)) |> 
  slice_head(20)

# Plot
vlcno_dc_hc <- volcano(tbl_dc_hc_bsln, top20_dc_hc_bsln) +
  annotation_custom(tableGrob(de_table(tbl_dc_hc_bsln, "DC vs. HC"),
                              theme = ttheme_minimal(
                                core = list(
                                  fg_params = list(hjust = 0.5))),
                              rows = NULL),
                    xmin = 1.6, xmax = 2.5,
                    ymin = 0.6, ymax = 3.5)

# Save plot
save_plot("output/figures/volcano_dc_vs_hc.pdf",
          vlcno_dc_hc,
          base_asp = 1.2, 
          ncol = 1,
          device = cairo_pdf)

## Figure 3b ####
# GSEA of pooled t-statistics from all disease groups vs. healthy controls,
# compare against the Reactome database
### Data wrangling ####
tbl_common_drh_vs_cntrl <- list("CDI vs. HC" = tbl_cdi_hc_bsln,
                                       "IVD vs. HC" = tbl_ibd_hc_bsln,
                                       "GDH+ vs. HC" = tbl_gdh_hc_bsln,
                                       "DC vs. HC" = tbl_dc_hc_bsln)

for (df in seq_along(tbl_common_drh_vs_cntrl)) {
  tbl_common_drh_vs_cntrl[[df]] <- tbl_common_drh_vs_cntrl[[df]] |>
    mutate(contrast = names(tbl_common_drh_vs_cntrl)[[df]]) |> 
    select(PROBEID, ENTREZID, t, contrast)
}

# Calculate mean t statistic
tbl_common_drh_vs_cntrl <- bind_rows(tbl_common_drh_vs_cntrl) |> 
  summarise(t = mean(t), .by = c(PROBEID, ENTREZID))

### Reactome GSEA ####
# Reactome gene sets downloaded on 2022-09-30, from
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
fgsea_rctm <- gmtPathways("raw/c2.cp.reactome.v2022.1.Hs.entrez.gmt")

# Function to extract and sort t-statistics for use with fgsea
t_sort <- function(df) {
  t_stats <- df$t
  
  # Name each t-statistic with Entrez ID
  names(t_stats) <- df$ENTREZID
  
  t_stats <- sort(t_stats, decreasing = TRUE)
  
  return(t_stats)
}

# fgsea for Reactome gene sets
set.seed(20240726)
common_drh_rctm_fgsea <- fgsea(pathways = fgsea_rctm,
                                      stats = t_sort(tbl_common_drh_vs_cntrl),
                                      eps = 0.0)

# Identify top 20 pathways with adjusted p-value < 0.05
common_drh_rctm_fgsea_de <- common_drh_rctm_fgsea |>
  filter(padj < 0.05) |> 
  slice_min(order_by = padj,
            n = 20) |> 
  pull(pathway)

# Plot summary table of significant pathways
tbl_fgsea_rctm_common_drh <- plotGseaTable(
  pathways = fgsea_rctm[common_drh_rctm_fgsea_de],
  stats = t_sort(tbl_common_drh_vs_cntrl),
  fgseaRes = common_drh_rctm_fgsea)

# Export summary table as pdf
ggsave("output/figures/fgsea_reactome_diarrhoea-common_hc_top-20.pdf",
  tbl_fgsea_rctm_common_drh,
  width = 190,
  units = "mm"
)

## Figure 4 ####
### CDI vs. IBD ####
# Import data
tbl_cdi_ibd_bsln <- fread("output/tables/de_stats_cdi_vs_ibd.csv")

# Identify top20 upregulated genes for point labelling
top20_cdi_ibd_bsln <- tbl_cdi_ibd_bsln |> 
  arrange(desc(logFC)) |> 
  slice_head(20)

# Plot
vlcno_cdi_ibd <- volcano(tbl_cdi_ibd_bsln, top20_cdi_ibd_bsln) +
  annotation_custom(tableGrob(de_table(tbl_cdi_ibd_bsln, "CDI vs. IBD"),
                              theme = ttheme_minimal(
                                core = list(
                                  fg_params = list(hjust = 0.5))),
                              rows = NULL),
                    xmin = 0.5, xmax = 1.8,
                    ymin = 11, ymax = 17)

# Save plot
save_plot("output/figures/volcano_cdi_vs_ibd.pdf",
          vlcno_cdi_ibd,
          base_asp = 1.2, 
          ncol = 1,
          device = cairo_pdf)

### CDI vs. GDH+ ####
# Import data
tbl_cdi_gdh_bsln <- fread("output/tables/de_stats_cdi_vs_gdh.csv")

# Identify top20 upregulated genes for point labelling
top20_cdi_gdh_bsln <- tbl_cdi_gdh_bsln |> 
  arrange(desc(logFC)) |> 
  slice_head(20)

# Plot
vlcno_cdi_gdh <- volcano(tbl_cdi_gdh_bsln, top20_cdi_gdh_bsln) +
  annotation_custom(tableGrob(de_table(tbl_cdi_gdh_bsln, "CDI vs. GDH"),
                              theme = ttheme_minimal(
                                core = list(
                                  fg_params = list(hjust = 0.5))),
                              rows = NULL),
                    xmin = 0.8, xmax = 1.2,
                    ymin = 8, ymax = 10)

# Save plot
save_plot("output/figures/volcano_cdi_vs_gdh.pdf",
          vlcno_cdi_gdh,
          base_asp = 1.2, 
          ncol = 1,
          device = cairo_pdf)

### CDI vs. DC ####
# Import data
tbl_cdi_dc_bsln <- fread("output/tables/de_stats_cdi_vs_dc.csv")

# Identify top20 upregulated genes for point labelling
top20_cdi_dc_bsln <- tbl_cdi_dc_bsln |> 
  arrange(desc(logFC)) |> 
  slice_head(20)

# Plot
vlcno_cdi_dc <- volcano(tbl_cdi_dc_bsln, top20_cdi_dc_bsln) +
  annotation_custom(tableGrob(de_table(tbl_cdi_dc_bsln, "CDI vs. DC"),
                              theme = ttheme_minimal(
                                core = list(
                                  fg_params = list(hjust = 0.5))),
                              rows = NULL),
                    xmin = 0.8, xmax = 1.2,
                    ymin = 6, ymax = 10)

# Save plot
save_plot("output/figures/volcano_cdi_vs_dc.pdf",
          vlcno_cdi_dc,
          base_asp = 1.2, 
          ncol = 1,
          device = cairo_pdf)


## Figure 6a ####
# Heatmap of z-scores from Ingenuity Pathway Analysis
# Import IPA z-scores
ipa_z <- fread("raw/Z-scores_CDIcomparison_canonical pathways.txt",
               na.string = ("N/A")) |>
  clean_names() |> 
  head(20)

# Save pathways in same order as IPA
ipa_levels <- ipa_z$canonical_pathways

colnames(ipa_z) <- sub("_adjp0.*", "", colnames(ipa_z))

# Wrangle for plotting
ipa_z_long <- ipa_z |> 
  pivot_longer(starts_with("ipa_"),
               names_to = "comparison",
               values_to = "z_score") |> 
  mutate(comparison = case_when(
    comparison == "ipa_cd_ivs_hc" ~ "CDI vs. HC",
    comparison == "ipa_cd_ivs_dc" ~ "CDI vs. DC",
    comparison == "ipa_cd_ivs_ibd" ~ "CDI vs. IBD",
    comparison == "ipa_cd_ivs_gdh" ~ "CDI vs. GDH",
    TRUE ~ NA_character_
  )) |> 
  mutate(comparison = factor(comparison,
                             levels = c(
                               "CDI vs. HC",
                               "CDI vs. DC",
                               "CDI vs. IBD",
                               "CDI vs. GDH"
                             ))) |>
  mutate(canonical_pathways = factor(canonical_pathways,
                                     levels = ipa_levels))

# Heatmap
heatmap_ipa_z <- ipa_z_long |> 
  ggplot(aes(x = comparison, y = canonical_pathways, fill = z_score)) +
  geom_tile(colour = "black") +
  scale_fill_gradientn(colors = rev(hcl.colors(255, "RdBu"))) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL,
       fill = "Z-score",
       y = NULL) +
  theme_minimal(base_size = 16) +
  theme(text = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,
                                   hjust = 0, vjust = 1))

# Export
save_plot("output/figures/heatmap_ipa_cdi_vs_drg_groups_z-scores.pdf",
          heatmap_ipa_z,
          ncol = 1.2,
          nrow = 1.5)

## Supplementary Figure 1 ####
# Incidence of C. diff. infection over time in various nations
### Wrangle data ####
cdi_inc <- data.frame(
  stringsAsFactors = FALSE,
  check.names = FALSE,
  Year = c("UK (per 100,000 persons)",
           "USA (per 100,000 persons)",
           "EU (per 1000 discharges or admissions)",
           "Australia (per 10,000 patient bed-days)"),
  `2011` = c(32.1, 140.9, NA, NA),
  `2012` = c(26.2, 145.8, NA, 4.3),
  `2013` = c(23.8, 141.8, NA, 3.94),
  `2014` = c(25.3, 141.6, NA, 3.81),
  `2015` = c(25.2, 148.6, NA, 3.85),
  `2016` = c(22.9, 146.3, 1.96, 3.91),
  `2017` = c(23.7, 130.3, 2.11, 3.92),
  `2018` = c(21.9, 130.1, 1.92, 3.97),
  `2019` = c(23.6, 121.2, 1.89, 4),
  `2020` = c(22.3, 101.3, 3.16, 3.95),
  `2021` = c(25.5, 110.2, NA, 4.82),
  `2022` = c(26.1, NA, NA, NA),
  `Jan-Apr.2023` = c(26L, NA, NA, NA)
)

# Pivot to long form for plotting
cdi_inc <- cdi_inc |>
  rename("location" = "Year") |> 
  pivot_longer(cols = where(is.numeric),
               names_to = "year",
               values_to = "incidence") |> 
  # Create numeric year column
  mutate(year_num = as.numeric(gsub("^Jan-Apr.", "", year))) |> 
  # Define units for plotting
  mutate(unit = case_when(
    grepl("10,000", location) ~ "Per 10,000 patient bed-days",
    grepl("100,000", location) ~ "Per 100,000 persons",
    TRUE ~ NA_character_
  )) |> 
  # Make copy of location with units of incidence
  mutate(location_detailed = location) |>
  mutate(location = gsub(
    " \\(per 100,000 persons\\)",
    "", location
  )) |> 
  mutate(location = gsub(
    " \\(per 10,000 patient-days\\)",
    "", location
  )) |> 
  mutate(location = gsub(
    " \\(per 10,000 patient bed-days\\)",
    "", location
  ))

### Plot ####
#### Incidence faceted by location
line_cdi_inc_facet_loc <- cdi_inc |> 
  filter(!is.na(incidence)) |> 
  mutate(label = if_else(year_num == max(year_num), location, NA_character_),
         .by = location) |>  
  ggplot(aes(x = year_num, y = incidence, colour = location)) +
  geom_point() +
  geom_line() +
  facet_wrap(~location_detailed, scales = "free_y",
             nrow = 1) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  scale_colour_manual(values = okabe_ito) +
  labs(x = NULL,
       y = "Incidence") +
  theme_bw_black(base_size = 18,
                 base_rect_size = 1,
                 base_family = "Helvetica") +
  theme(legend.position = "none") +
  background_grid(major = "y")

# Save plot
save_plot("output/figures/cdi_incidence_faceted_by_location.pdf",
          line_cdi_inc_facet_loc,
          base_asp = 1.2,
          ncol = 4)

## Supplementary Figure 2 ####
# Bar chart of participants with each condition, split among sites
site_data <- fread("raw/deidentified_site_by_condition.csv") |> 
  mutate(condition_labels = factor(condition_labels,
                                   levels = c("CDI", "GDH", "IBD",  "DC", "HC")))

bar_site_condition <- site_data |>
  summarise(count = n(),
            .by = c(condition_labels, site_label)) |>
  ggplot(aes(x = site_label, y = count, fill = condition_labels)) +
  geom_col(position = position_dodge(preserve = "single"), # Consistent width
           colour = "black") +
  scale_fill_manual(values = okabe_ito) +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  theme_bw_black(base_size = 14) +
  background_grid(major = "y", minor = "y") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.93, 0.74),
    legend.box.background = element_rect(colour = "black", linewidth = 1),
    legend.text = element_text(size = 9)
  ) +
  labs(x = NULL, y = "Count", fill = NULL)

# Save plot
save_plot("output/figures/site_by_condition_bar_chart.svg",
          bar_site_condition,
          base_asp = 2,
          ncol = 1)

## Supplementary Figure 3 ####
# Format data
participants_setting <- data.frame(
  stringsAsFactors = FALSE,
  check.names = FALSE,
  Group = c("CDI", "GDH", "IBD", "DC", "HC"),
  `Hospital ward` = c(77L, 34L, 37L, 45L, NA),
  `Critical Care Unit` = c(1L, 2L, NA, NA, NA),
  `Day Case Unit` = c(NA, 1L, 1L, NA, 12L),
  `Outpatient clinics` = c(NA, NA, 2L, NA, 31L),
  `Other` = c(NA, NA, NA, NA, 8L)
) |> 
  pivot_longer(cols = where(is.numeric),
               names_to = "setting",
               values_to = "count")

participants_setting$Group <- factor(participants_setting$Group,
       levels = c("CDI", "GDH", "IBD",  "DC", "HC"))

# Bar chart of patients recruited per setting
bar_particpants_per_setting <- participants_setting |> 
  ggplot(
    aes(x = Group,
        y = count,
        fill = setting
    )
  ) +
  geom_col(
    position = position_dodge(preserve = "single"), # Consistent width
    colour = "black"
  ) +
  scale_fill_manual(
    values = c("#134074", "#bfab25", "#4ea699", "#efb0a1", "#df2935")) +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  theme_bw_black(base_size = 14) +
  background_grid(major = "y", minor = "y") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.882, 0.77),
        legend.box.background = element_rect(colour = "black", linewidth = 1),
        legend.text = element_text(size = 9)
  ) +
  labs(x = NULL,
       y = "Count",
       fill = NULL)

# Save plot
save_plot("output/figures/particpants_per_setting_bar_chart.svg",
          bar_particpants_per_setting,
          base_asp = 2,
          ncol = 1)

## Supplementary Figure 4 ####
# Bar chart of ribotyped stool samples
### Wrangle data ####
ribotype <- data.frame(
  stringsAsFactors = FALSE,
  ribotype = c("001",
               "002","005","012","013","014","015",
               "017","018","020","023","026","027","046",
               "049","050","078","087","Other",
               "Non-toxigenic strain","Negative culture",
               "No sample"),
  CDI = c(2L,9L,
          3L,2L,NA,9L,2L,2L,NA,4L,9L,NA,4L,
          NA,1L,NA,9L,1L,11L,NA,4L,6L),
  GDH = c(NA,2L,
          NA,NA,1L,3L,1L,NA,1L,1L,NA,2L,2L,
          2L,NA,1L,2L,NA,2L,7L,3L,8L)
) |> 
  pivot_longer(cols = c(CDI, GDH),
               names_to = "cohort",
               values_to = "count")

### Plot ####
# Bar chart
bar_ribotype <- ribotype |> 
  ggplot(aes(x = ribotype, y = count, fill = cohort)) +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_manual(values = okabe_ito[c(1, 3)]) +
  scale_y_continuous(expand = expansion(c(0, 0.1)),
                     breaks = scales::breaks_pretty()) +
  theme_bw_black(base_size = 12) +
  background_grid(major = "y", minor = "y") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.05, 0.91),
        legend.background =  element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(10, "pt"),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1, vjust = 1)) +
  labs(x = NULL, y = "Count",
       fill = NULL)

# Save plot
save_plot("output/figures/ribotype_bar_chart.pdf",
          bar_ribotype,
          base_asp = 2,
          ncol = 1)

## Supplementary Figure 5 ####
# Formate data
gdh_neg_cdt_neg_ribotype <- data.frame(
  stringsAsFactors = FALSE,
  Ribotype = c("001", "005", "081", "106", "Other"),
  DC = c(3L, 1L, NA, 1L, 2L),
  IBD = c(NA, NA, 1L, NA, 1L),
  HC = c(1L, NA, NA, NA, NA)
) |> 
  pivot_longer(cols = c(DC,, IBD, HC),
               names_to = "cohort",
               values_to = "count") |> 
  mutate(cohort = factor(cohort, levels = c("IBD",  "DC", "HC")))

### Plot ####
# Bar chart
bar_gdh_neg_cdt_neg_ribotype <- gdh_neg_cdt_neg_ribotype |> 
  ggplot(aes(x = Ribotype, y = count, fill = cohort)) +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_manual(values = okabe_ito[c(3, 4, 5)]) +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  theme_bw_black(base_size = 16) +
  background_grid(major = "y", minor = "y") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.95, 0.91),
        legend.background =  element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(10, "pt")) +
  labs(x = NULL, y = "Count",
       fill = NULL)

# Save plot
save_plot("output/figures/gdh_neg_cdt_neg_ribotype_bar_chart.pdf",
          bar_gdh_neg_cdt_neg_ribotype,
          base_asp = 2,
          ncol = 1)

## Supplementary Figure 6 ####
### Import data ####
# Import prcomp object from RMA-normalised data filtered to Entrez genes
pca_rma_genes <- readRDS("processed/pca_rma_entrez.RDS")

# Microarray expressionSet
eset <- readRDS("processed/RMAExpressionSet.RDS")
eset$condition <- factor(eset$condition,
                         levels = c("CDI", "GDH", "IBD", "DC", "HC"))


# Create new condition labels and reorder
eset$condition_labels <- if_else(eset$condition == "DC GDH+",
                                 "GDH", eset$condition)
eset$condition_labels <- factor(eset$condition_labels,
                             levels = c("CDI", "GDH", "IBD",  "DC", "HC"))
### Plotting functions ####
# Function to return PC values for plotting
pca_df <-
  function(pca_obj, metadata, group, pc = c(1L, 2L),
           pca_object = FALSE,
           large_palette = FALSE, id = FALSE,
           discrete = TRUE) {
    
    # Check all required inputs are present
    if (missing(pca_obj)) {
      stop("PCA object must be supplied")
    }
    if (missing(metadata)) {
      stop("Metadata must be supplied")
    }
    if (missing(group)) {
      stop("Group must be defined")
    }
    
    # Function to detect whole numbers, from examples in ?is.integer
    is.wholenumber <-
      function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
    # Check that PCs to plot are integers
    if (any(!is.wholenumber(as.numeric((pc))))) {
      stop("PCs to be plotted must be numeric integers")
    }
    
    # Define objects and metadata
    metadat <- Biobase::pData(metadata)
        grp <- metadat[, group]
    
        
    # Define principal components to be plotted
    pca_raw <- pca_obj
    selected_pc1 <- as.integer(pc[1])
    selected_pc2 <- as.integer(pc[2])
    
    if (pca_object == FALSE) {
    pc1 <-  pca_raw$x[, selected_pc1]
    pc2 <-  pca_raw$x[, selected_pc2]
    } else {
      # For objects of class 'PCA'
      pc1 <- pca_raw[["rotated"]][, selected_pc1]
      pc2 <- pca_raw[["rotated"]][, selected_pc2]
    }
    
    # If ID columns to add to exported data are required
    if (any(id != FALSE)) {
      id_label <- metadat[, id]
    } else {
      id_label <- NA_character_
    }
    
    # Prepare data.frame to plot
    data_gg <- data.frame(pc1 = pc1,
                          pc2 = pc2,
                          grp = grp)
    # Add ID columns
    data_gg <- cbind.data.frame(data_gg, id_label)
    
    # Extract PCA values for plot axis labels
    percent_var <-
      round(100 * pca_raw$sdev ^ 2 / sum(pca_raw$sdev ^ 2), 1)
    
    data_gg <- cbind.data.frame(data_gg,
                                "selected_pc1_var" = percent_var[selected_pc1],
                                "selected_pc2_var" = percent_var[selected_pc2])
    
    return(data_gg)
  }

# Format data for PC scatterplot
pca_dat_genes_cond <- pca_df(pca_rma_genes, eset, "condition_labels",
                             pca_object = FALSE)

# Scatterplot of PC1 vs. PC2 from gene-level data, coloured by condition
scttr_pca_all_genes <- pca_dat_genes_cond |>
  ggplot(aes(pc1, pc2)) +
  geom_point(aes(fill = grp),
             shape = 21, colour = "black",
             alpha = 0.8) +
  stat_ellipse(aes(colour = grp),
               alpha = 0.5) +
  xlab(paste0("PC1 (Variance explained: ",
              max(pca_dat_genes_cond$selected_pc1_var), "%)")) +
  ylab(paste0("PC2 (Variance explained: ",
              max(pca_dat_genes_cond$selected_pc2_var), "%)")) +
  scale_fill_manual(values = okabe_ito) +
  scale_colour_manual(values = okabe_ito) +
  theme_bw_black(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  labs(fill = NULL) +
  guides(colour = "none")

# Save plot
save_plot("output/figures/pca_scatterplot_by_condition.svg",
          scttr_pca_all_genes,
          base_asp = 1,
          ncol = 1)

## SessionInfo ####
writeLines(capture.output(sessionInfo()),
           "output/sessionInfo/03_manuscript_figures_sessionInfo.txt")

