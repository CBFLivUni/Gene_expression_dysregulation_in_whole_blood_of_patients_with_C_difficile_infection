## Packages used ####
library(tidytable)
library(oligo)
library(stringr)
library(clariomdhumantranscriptcluster.db)

## Import data ####
### Metadata ####
meta <- fread("raw/metadata.csv")
colnames(meta) <- gsub("characteristics: ", "", colnames(meta))

# Reformat CEL file names to match GEO
meta$`CEL file` <- gsub("_\\(|\\)", "_", meta$`CEL file`)

### Microarray data ####
# Specify file path to CEL files obtained from GEO GSE276395
cel_files <- list.celfiles("raw/GSE276395_RAW",
                           full.names = TRUE,
                           listGzipped = TRUE)

# Check that file order matches metadata
cel_files_trunc <- str_extract(cel_files,
                               "(?<=raw/GSE276395_RAW/GSM\\d{7}_).+(?=.gz)")
# Check order matches, should return TRUE
identical(cel_files_trunc, meta$`CEL file`)
rm(cel_files_trunc)

# Import microarray data as ExpressionSet object
eset_raw <- read.celfiles(cel_files,
                      verbose = FALSE)
rm(cel_files)

# Add metadata
meta <- as.data.frame(
  meta,
  row.names = list.celfiles("raw/GSE276395_RAW", listGzipped = TRUE))
phenoData(eset_raw) <- AnnotatedDataFrame(meta)

## RMA processing ####
# Note that if running on an HPC, to avoid an error run the following in R:
# BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
eset <- oligo::rma(eset_raw, target = "core")

## Filter transcription clusters by intensity ####
intnsty_cutoff <- 3

# Transcript clusters (rows) are filtered based on at least 37 samples having
# expression greater than the threshold.
smpl_num_cutoff <- 37

# Index of transcript clusters meeting the above criteria
idx <- rowSums(Biobase::exprs(eset) > intnsty_cutoff) >= smpl_num_cutoff

# Check number of rows included and excluded after this step
table(idx)

# Filter data
eset <- subset(eset, idx)

# Check number of rows matches expected number of rows
# Should return TRUE
nrow(eset) == table(idx)[[2]]
rm(smpl_num_cutoff, idx, intnsty_cutoff)

## Annotate transcription clusters ####
### Match transcript clusters (TCs) to Entrez gene IDs ####
# Create data frame with Entrez gene IDs from clariomdhumantranscriptcluster.db
eset_anno <- AnnotationDbi::select(clariomdhumantranscriptcluster.db,
                                  keys = featureNames(eset),
                                  keytype = "PROBEID",
                                  columns = c("SYMBOL", "GENENAME", "ENTREZID"))

# Remove rows with unassigned gene symbols
eset_anno <- filter(eset_anno, !is.na(SYMBOL))

# Identify transcript clusters assigned to multiple genes
multi_matches <- eset_anno |> 
  summarise(no_of_matches = n_distinct(SYMBOL), .by = PROBEID) |> 
  filter(no_of_matches > 1)
nrow(multi_matches) # 1049 with multiple genes mapped

# Create indices of probes to remove from expression data based on multi_matches
idx_exclude <- (featureNames(eset) %in% multi_matches$PROBEID)
rm(multi_matches)

### Annotate data with gene symbols ####
# Remove TC IDs with multiple gene matches from data
eset <- subset(eset, !idx_exclude)
validObject(eset)
rm(idx_exclude)

# Add column of PROBEIDs to feature data
fData(eset)$PROBEID <- rownames(fData(eset))

# Match gene ID annotations to feature data
fData(eset) <- dplyr::left_join(fData(eset), eset_anno, by = "PROBEID")

# Restore rownames
rownames(fData(eset)) <- fData(eset)$PROBEID
validObject(eset)
rm(eset_anno)

### Filter transcript clusters with no Entrez match  ####
eset <- subset(eset, !is.na(eset@featureData@data$SYMBOL))

### Filter transcript probes with multiple Entrez matches ####
# Within each gene ID, select the TC with the greatest mean expression

# Calculate mean expression for each TC
tc_mean <- rowMeans(exprs(eset))

# Obtain order index of expression, from greatest to smallest
tc_ordered <- order(tc_mean, decreasing = TRUE)

# Reorder by TC mean expression, highest to lowest
eset <- eset[tc_ordered, ]
rm(tc_mean, tc_ordered)

# Identify duplicated gene IDs to remove
# Due to ordering by mean expression, on average these TCs will have lower
# expression
tc_to_remove <- duplicated(fData(eset)$ENTREZID)

# Remove these TCs
eset <- eset[!tc_to_remove, ]
rm(tc_to_remove)

## Save filtered ExpressionSet object ####
saveRDS(eset, "processed/RMAExpressionSet.RDS")

## SessionInfo ####
writeLines(capture.output(sessionInfo()),
           "output/sessionInfo/01_process_cel_sessionInfo.txt")


