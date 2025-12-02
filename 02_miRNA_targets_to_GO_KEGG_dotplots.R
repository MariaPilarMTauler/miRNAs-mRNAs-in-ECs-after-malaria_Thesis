###############################################
# From miRNA targets (multiMiR) to GO / KEGG enrichment
#
# - Input: table of miRNA → target gene (symbols)
# - Output: GO / KEGG dotplots for targets of upregulated (or selected) miRNAs
###############################################

# Load packages #
library(readxl)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

## 1) Load multiMiR target table #
## This can be the same "targets_df" you export in the Sankey script: ##
##   mature_miRNA_id, target_symbol, etc.

targets_file <- "Sankey_Master_Table_targets_only.xlsx"  # example
targets_sheet <- "targets_df"                            # or use CSV

## Example if you exported as Excel with a sheet "targets_df" ##
targets_df <- readxl::read_excel(targets_file, sheet = targets_sheet)

### Make sure we have: ###
#   - mature_miRNA_id
#   - target_symbol
targets_df <- targets_df %>%
  dplyr::rename(
    miRNA  = mature_miRNA_id,
    SYMBOL = target_symbol
  )

##### 2) (Optional) restrict to targets of specific miRNAs ###
## Example: only targets of upregulated miRNAs of interest ##
## You can skip this if you want enrichment on all targets ##

## Vector of miRNAs of interest (e.g. upregulated) ##
## mirna_of_interest <- c("hsa-miR-193a-3p", "hsa-miR-485-5p")  ##

targets_sel <- targets_df %>%
  filter(miRNA %in% mirna_of_interest)

## Get unique target gene symbols ##
target_genes <- unique(targets_sel$SYMBOL)

cat("Number of unique target genes:", length(target_genes), "\n")

## 3) Map target gene symbols to Entrez IDs ##
gene2entrez <- bitr(
  target_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

entrez_targets <- unique(gene2entrez$ENTREZID)

## 4) GO Biological Process enrichment ##
ego_targets <- enrichGO(
  gene          = entrez_targets,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

## 5) KEGG enrichment ##
ekegg_targets <- enrichKEGG(
  gene          = entrez_targets,
  organism      = "hsa",
  keyType       = "kegg",
  pvalueCutoff  = 0.05
)

ekegg_targets <- setReadable(ekegg_targets, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# 6) Dotplots for miRNA target enrichment #
out_dir <- "GO_KEGG_dotplots_miRNA_targets"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# GO BP dotplot for miRNA targets
p_go_targets <- dotplot(ego_targets, showCategory = 20) +
  ggtitle("GO BP enrichment - miRNA target genes")

ggsave(
  filename = file.path(out_dir, "GO_BP_miRNA_targets_dotplot.png"),
  plot     = p_go_targets,
  width    = 7,
  height   = 5,
  dpi      = 300
)

### KEGG dotplot for miRNA targets ##
p_kegg_targets <- dotplot(ekegg_targets, showCategory = 20) +
  ggtitle("KEGG enrichment - miRNA target genes")

ggsave(
  filename = file.path(out_dir, "KEGG_miRNA_targets_dotplot.png"),
  plot     = p_kegg_targets,
  width    = 7,
  height   = 5,
  dpi      = 300
)

cat("✅ GO / KEGG dotplots for miRNA targets saved in:", out_dir, "\n")
