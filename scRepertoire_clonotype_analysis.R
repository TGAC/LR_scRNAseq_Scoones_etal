# ===========================================================
# Clonal Analysis of PBMC Samples Using scRepertoire
# ===========================================================
# Description:
# This script processes and analyses TCR/BCR repertoires from single-cell RNA sequencing (scRNA-seq) data 
# across multiple long-read sequencing platforms using scRepertoire.
#
# Requirements:
# - This script assumes that users have:
#   1. Downloaded the Seurat objects (*.rds files) for each sequencing platform.
#   2. Obtained TRUST4 barcode report files (*.tsv) for each platform.
#   3. Installed the required R packages: ggplot2, dplyr, scRepertoire.
#
# Outputs:
# - Clonal analysis metrics including overlap, homeostasis, diversity, and abundance comparisons.
# - Figures visualising TCR and BCR clonal distributions across sequencing platforms.
# ===========================================================

rm(list=ls())

#Required Libraries
library(ggplot2)
library(dplyr)
library(scRepertoire)

# ===========================================================
# Define Input File Paths
# ===========================================================

seurat_files <- list(
  ONT = "2025_ONT_1b2a_sobj_publication.rds",
  revio = "2025_PacBio_MASseq_Revio_1b2a_sobj_publication.rds",
  revio_jc = "2025_PacBio_MASseq_Revio_Jumpcode_1b2a_sobj_publication.rds",
  sequel = "2025_PacBio_MASseq_SequelIIe_1b2a_sobj_publication.rds",
  sequel_jc = "2025_PacBio_MASseq_SequelIIe_Jumpcode_1b2a_sobj_publication.rds"
)

trust4_files <- list(
  ONT = "path/to/nanopore_TRUST4/sample1B_ONT_barcode_report.tsv",
  revio = "path/to/revio_TRUST4/sample1B/sample1B_revio_barcode_report.tsv",
  revio_jc = "path/to/revio_jumpcode_TRUST4/sample1B/sample1B_reviojumpcode_barcode_report.tsv",
  sequel = "path/to/sequel_TRUST4/sample1B/sample1B_sequel_barcode_report.tsv",
  sequel_jc = "path/to/sequel_jumpcode_TRUST4/sample1B/sample1B_sequeljumpcode_barcode_report.tsv"
)

# ===========================================================
# Load Seurat Objects and Extract Barcodes
# ===========================================================

seurat_objects <- lapply(seurat_files, readRDS)

# Function to clean barcodes from Seurat objects
clean_seurat_barcodes <- function(barcodes) {
  barcodes <- gsub("-1_1$", "", barcodes)          # Remove "-1_1" suffix (PacBio Revio & Revio JC)
  barcodes <- gsub("-1_MAS_1$", "", barcodes)      # Remove "-1_MAS_1" suffix (PacBio Sequel)
  barcodes <- gsub("-1_MAS_JC_1$", "", barcodes)   # Remove "-1_MAS_JC_1" suffix (PacBio Sequel JC)
  barcodes <- sub("-1$", "", barcodes)            # Remove simple "-1" suffix (ONT, 10X)
  return(barcodes)
}

# Extract and clean barcodes
valid_barcodes <- lapply(seurat_objects, function(obj) clean_seurat_barcodes(colnames(obj)))

# ===========================================================
# Load TRUST4 Barcode Reports
# ===========================================================

trust4_data <- lapply(trust4_files, function(file) read.csv(file, sep = '\t'))

# Function to clean TRUST4 barcodes
clean_trust4_barcodes <- function(df) {
  if (!"X.barcode" %in% colnames(df)) stop("Barcode column missing in TRUST4 data")
  
  colnames(df)[colnames(df) == "X.barcode"] <- "barcode"
  df$barcode <- sub("^.*_", "", df$barcode)  # Remove sample prefixes
  df$barcode <- gsub("-1$", "", df$barcode)  # Remove "-1" suffix if present
  return(df)
}

# Apply barcode cleaning to TRUST4 data
trust4_data_cleaned <- lapply(trust4_data, clean_trust4_barcodes)

# ===========================================================
# Filter TRUST4 Data to Retain Only Seurat-Matching Barcodes
# ===========================================================

trust4_filtered <- mapply(function(df, sample) {
  df <- df[df$barcode %in% valid_barcodes[[sample]], ]
  return(df)
}, trust4_data_cleaned, names(valid_barcodes), SIMPLIFY = FALSE)

# ===========================================================
# Load TRUST4 Data into scRepertoire
# ===========================================================

T4data <- loadContigs(input = trust4_filtered, format = "TRUST4")
#head(T4data)

# ===========================================================
# TCR and BCR Clonal Analysis
# ===========================================================

# Filter for TCR chains (TRA, TRB, TRG, TRD)
T4data_TCR <- lapply(T4data, function(df) {
  df %>% filter(chain %in% c("TRA", "TRB", "TRG", "TRD"))  # Keep only TCR chains
})

# Run combineTCR with filtered TCR data
combined.TCR <- combineTCR(T4data_TCR, 
                           samples = c("ONT", "sequel", "sequelJC", "revio", "revioJC"), 
                           removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)

head(combined.TCR)

# Filter for BCR chains (IGH, IGK, IGL)
T4data_BCR <- lapply(T4data, function(df) {
  df %>% filter(chain %in% c("IGH", "IGK", "IGL"))
})

# Run combineBCR with filtered TCR data
combined.BCR <- combineBCR(T4data_BCR, samples = names(valid_barcodes), threshold = 0.85)


head(combined.BCR)

# ===========================================================
# Clonal Visualisation
# ===========================================================
TCRclones <- clonalCompare(combined.TCR, 
                           top.clones = 10, 
                           samples = c("ONT", "sequel", "sequelJC", "revio", "revioJC"), 
                           cloneCall="strict", 
                           graph = "alluvial", 
                           relabel.clones = TRUE)

BCRclones <- clonalCompare(combined.BCR, 
                           top.clones = 10, 
                           samples = c("ONT", "sequel", "sequelJC", "revio", "revioJC"), 
                           cloneCall="strict", 
                           graph = "alluvial", 
                           relabel.clones = TRUE)

TCRclones | BCRclones

# ===========================================================
# Subset Analysis for αβT and γδT Cells
# ===========================================================

# Function to filter data by cell type
filter_dataset <- function(dataset, celltype) {
  if (!"cell_type" %in% colnames(dataset)) {
    stop("Error: 'cell_type' column not found in dataset")
  }
  dataset %>% filter(cell_type == celltype)
}

# Generate subsets
T4list_abT <- lapply(T4data, filter_dataset, celltype = "abT")
T4list_gdT <- lapply(T4data, filter_dataset, celltype = "gdT")

# Process αβT subset
T4data_abT <- loadContigs(input = T4list_abT, format = "TRUST4")

# May need to filter to remove invalid "Non" chain values
T4data_TCR_abT <- lapply(T4data_abT, function(df) {
  df %>% filter(chain %in% c("TRA", "TRB", "TRG", "TRD"))  # Keep only valid TCR chains
})

# Run combineTCR
combined.TCR.abT <- combineTCR(T4data_TCR_abT, 
                               samples = c("ONT", "sequel", "sequelJC", "revio", "revioJC"), 
                               removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)

# Visualise abT clones
abTclones <- clonalCompare(combined.TCR.abT, 
                           top.clones = 10, 
                           samples = c("ONT", "sequel", "sequelJC", "revio", "revioJC"), 
                           cloneCall="strict", 
                           graph = "alluvial", 
                           relabel.clones = TRUE)
abTclones

# Process γδT subset
T4data_gdT <- loadContigs(input = T4list_gdT, format = 'TRUST4')

# May need to filter to remove invalid "Non" chain values
T4data_TCR_gdT <- lapply(T4data_gdT, function(df) {
  df %>% filter(chain %in% c("TRA", "TRB", "TRG", "TRD"))  # Keep only valid TCR chains
})

# Run combineTCR for gdT cells
combined.TCR.gdT <- combineTCR(T4data_TCR_gdT, 
                               samples = c("ONT", "sequel", "sequelJC", "revio", "revioJC"), 
                               removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)

# Visualise gdT clones
gdTclones <- clonalCompare(combined.TCR.gdT, 
                           top.clones = 10, 
                           samples = c("ONT", "sequel", "sequelJC", "revio", "revioJC"), 
                           cloneCall="strict", 
                           graph = "alluvial", 
                           relabel.clones = TRUE)

abTclones | gdTclones | BCRclones

# ===========================================================
# Visualisation of clonal structure 
# ===========================================================

# Clonal scatter per cell type
clonalScatter(combined.TCR.abT, cloneCall ="gene", x.axis = "ONT", y.axis = "revio", dot.size = "total", graph = "proportion")
clonalScatter(combined.TCR.gdT, cloneCall ="gene", x.axis = "ONT", y.axis = "revio", dot.size = "total", graph = "proportion")
clonalScatter(combined.BCR, cloneCall ="gene", x.axis = "ONT", y.axis = "revio", dot.size = "total", graph = "proportion")

# Total numbers of unique clones using “both” for combined chain visualization (can modify  “TRA”, “TRB”, “TRD”, “TRG”, “IGH” or “IGL” to select single chain)

clonalQuant(combined.BCR,     cloneCall="strict",    chain = "both",     scale = FALSE)
clonalQuant(combined.TCR.abT,  cloneCall="strict",  chain = "both", scale = FALSE)
clonalQuant(combined.TCR.gdT,   cloneCall="strict",  chain = "both",  scale = FALSE)

#  Using the the abundance of clones across groupings estimate the rarefaction across groupings. 
#Underlying the rarefaction calculation is the use of observed receptor of abundance to compute diversity where hill.number = 1, Shannon Diversity (q = 1), and plot.type = 2 - sample completeness curve

#BCR 
clonalRarefaction(combined.BCR,
                  plot.type = 2,
                  hill.numbers = 1,
                  n.boots = 2)

#gdT 
clonalRarefaction(combined.TCR.gdT,
                  plot.type = 2,
                  hill.numbers = 1,
                  n.boots = 2)
#abT
clonalRarefaction(combined.TCR.abT,
                  plot.type = 2,
                  hill.numbers = 1,
                  n.boots = 2)

# ===========================================================
# Clonal Homeostasis
# ===========================================================
#Relative space occupied by clones at specific proportions ie. what percentage of each run is filled by clones in distinct proportions - proportional cut points for cloneSize used are c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded = 1).
clonalHomeostasis(combined.BCR,   cloneCall = "gene",  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded = +  1))
clonalHomeostasis(combined.TCR.abT,  cloneCall = "gene", cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded = + 1))
clonalHomeostasis(combined.TCR.gdT, cloneCall = "gene",  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded = 1))

# ===========================================================
# Clonal Overlap Analysis
# ===========================================================

#raw n of overlapping clones
clonalOverlap(combined.TCR.abT, cloneCall = "strict", method = "raw")
clonalOverlap(combined.TCR.gdT, cloneCall = "strict", method = "raw")
clonalOverlap(combined.BCR, cloneCall = "strict", method = "raw")

# weighted overlap metric that accounts for both shared and abundant clonotypes - quantify how similarly two platforms capture clonal diversity
clonalOverlap(combined.TCR.gdT, cloneCall = "strict", method = "morisita")
clonalOverlap(combined.BCR, cloneCall = "strict", method = "morisita")
clonalOverlap(combined.TCR.abT, cloneCall = "strict", method = "morisita")

# ===========================================================
# End of Script
# ===========================================================
