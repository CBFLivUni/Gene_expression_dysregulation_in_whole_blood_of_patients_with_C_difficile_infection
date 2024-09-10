## Set up ####
library(oligo)
library(limma)
library(splines)
library(tidytable)

## Import data ####
# ExpressionSet object
eset <- readRDS("processed/RMAExpressionSet.RDS")

## Principal Componenent Analysis ####
# PCA of RMA data with annotated genes ####
expr_eset_genes <- Biobase::exprs(eset)
pca_genes <- prcomp(t(expr_eset_genes), scale = FALSE)

# Export prcomp object
saveRDS(pca_genes,
        "processed/pca_rma_entrez.RDS")

## Differential expression analysis ####
### Diarrhoea groups vs. healthy controls ####
# Filter expression data to baseline timepoint samples only
eset_bsln <- eset[, eset$time_point_sample_weeks == 0L]

# Check that all participant IDs are unique, should return TRUE
length(unique(eset_bsln$participant_id)) == ncol(eset_bsln)

#### Define model variables ####
# Define age
age <- ns(eset_bsln$age, df = 3)

# Define sex
sex <- factor(eset_bsln$sex)

# Define disease conditions
condition <- eset_bsln$condition

#### CDI vs. healthy controls ####
# Design matrix
design <- model.matrix(~ 0 + condition + age + sex)

# Rename design columns
colnames(design) <- c(
  "CDI", "DC", "GDH", "HC", "IBD", 
  "age1", "age2", "age3",
  "Male"
)

# Fit models
fit <- lmFit(eset_bsln, design = design)

# Define contrasts and identify DE genes
cntrsts <- makeContrasts(
  "CDI - HC",
  levels = design)

# Estimate coefficients and SEs for contrasts
fit_contrsts <- contrasts.fit(fit, contrasts = cntrsts)

# Compute moderated statistics
fit_contrsts <- eBayes(fit_contrsts)

# Identify DE genes for CDI vs. HC
tbl_cdi_hc <- topTable(fit_contrsts,
                       coef = "CDI - HC",
                       number = Inf) |> 
  mutate(pval_under_05 = ifelse(adj.P.Val < 0.05, TRUE, FALSE))

# Export csv file of CDI vs. healthy controls
fwrite(tbl_cdi_hc, "output/tables/de_stats_cdi_vs_hc.csv")

#### IBD vs. healthy controls ####
# Specify contrasts
cntrsts <- makeContrasts(
  "IBD - HC",
  levels = design)

# Estimate coefficients and SEs for contrasts
fit_contrsts <- contrasts.fit(fit, contrasts = cntrsts)

# Compute moderated statistics
fit_contrsts <- eBayes(fit_contrsts)

# Identify DE genes for CDI vs. healthy controls
tbl_ibd_hc <- topTable(fit_contrsts,
                       coef = "IBD - HC",
                       number = Inf) |> 
  mutate(pval_under_05 = ifelse(adj.P.Val < 0.05, TRUE, FALSE))

# Export csv file of IBD vs. healthy controls
fwrite(tbl_ibd_hc, "output/tables/de_stats_ibd_vs_hc.csv")

#### GDH+ diarrhoea vs. healthy controls ####
# Define contrasts and identify DE genes
cntrsts <- makeContrasts(
  "GDH - HC",
  levels = design)

# Estimate coefficients and SEs for contrasts
fit_contrsts <- contrasts.fit(fit, contrasts = cntrsts)

# Compute moderated statistics
fit_contrsts <- eBayes(fit_contrsts)

# Identify DE genes for CDI vs. healthy controls
tbl_gdh_hc <- topTable(fit_contrsts,
                       coef = "GDH - HC",
                       number = Inf) |> 
  mutate(pval_under_05 = ifelse(adj.P.Val < 0.05, TRUE, FALSE))

# Export csv file of GDH+ diarrhoea vs. healthy controls
fwrite(tbl_gdh_hc, "output/tables/de_stats_gdh_vs_hc.csv")

#### Diarrhoea control vs. healthy controls ####
# Define contrasts and identify DE genes
# Specify contrasts
cntrsts <- makeContrasts(
  "DC - HC",
  levels = design)

# Estimate coefficients and SEs for contrasts
fit_contrsts <- contrasts.fit(fit, contrasts = cntrsts)

# Compute moderated statistics
fit_contrsts <- eBayes(fit_contrsts)

# Identify DE genes for CDI vs. diarrhoea from other causes
tbl_dc_hc <- topTable(fit_contrsts,
                      coef = "DC - HC",
                      number = Inf) |> 
  mutate(pval_under_05 = ifelse(adj.P.Val < 0.05, TRUE, FALSE))

# Export csv file of diarrhoea controls vs. healthy controls
fwrite(tbl_dc_hc, "output/tables/de_stats_dc_vs_hc.csv")

### CDI vs. diarrhoea groups ####
##### CDI vs. IBD ####
# Define contrasts and identify DE genes
cntrsts <- makeContrasts(
  "CDI - IBD",
  levels = design)

# Estimate coefficients and SEs for contrasts
fit_contrsts <- contrasts.fit(fit, contrasts = cntrsts)

# Compute moderated statistics
fit_contrsts <- eBayes(fit_contrsts)

# Identify DE genes for CDI vs. IBD
tbl_cdi_ibd <- topTable(fit_contrsts,
                        coef = "CDI - IBD",
                        number = Inf) |> 
  mutate(pval_under_05 = ifelse(adj.P.Val < 0.05, TRUE, FALSE))

# Export csv file of CDI vs. IBD DE stats
fwrite(tbl_cdi_ibd, "output/tables/de_stats_cdi_vs_ibd.csv")

##### CDI vs. GDH ####
# Define contrasts and identify DE genes
cntrsts <- makeContrasts(
  "CDI - GDH",
  levels = design)

# Estimate coefficients and SEs for contrasts
fit_contrsts <- contrasts.fit(fit, contrasts = cntrsts)

# Compute moderated statistics
fit_contrsts <- eBayes(fit_contrsts)

# Identify DE genes for CDI vs. IBD
tbl_cdi_gdh <- topTable(fit_contrsts,
                        coef = "CDI - GDH",
                        number = Inf) |> 
  mutate(pval_under_05 = ifelse(adj.P.Val < 0.05, TRUE, FALSE))

# Export csv file of CDI vs. GDH DE stats
fwrite(tbl_cdi_gdh, "output/tables/de_stats_cdi_vs_gdh.csv")

##### CDI vs. DC ####
# Define contrasts and identify DE genes
cntrsts <- makeContrasts(
  "CDI - DC",
  levels = design)

# Estimate coefficients and SEs for contrasts
fit_contrsts <- contrasts.fit(fit, contrasts = cntrsts)

# Compute moderated statistics
fit_contrsts <- eBayes(fit_contrsts)

# Identify DE genes for CDI vs. DC
tbl_cdi_dc <- topTable(fit_contrsts,
                        coef = "CDI - DC",
                        number = Inf) |> 
  mutate(pval_under_05 = ifelse(adj.P.Val < 0.05, TRUE, FALSE))

# Export csv file of CDI vs. IBD DE stats
fwrite(tbl_cdi_dc, "output/tables/de_stats_cdi_vs_dc.csv")

## SessionInfo ####
writeLines(capture.output(sessionInfo()),
           "output/sessionInfo/02_pca_dea_sessionInfo.txt")

