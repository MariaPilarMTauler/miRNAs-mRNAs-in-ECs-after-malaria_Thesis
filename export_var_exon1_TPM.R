library(dplyr)
library(readr)
library(writexl)

# === INPUT FILES ===
tpm_file <- "counts_tpm.txt"
var_gene_file <- "HB3_var_gene_list.txt"  # << updated to new clean list
output_file <- "var_genes_exon1_TPM.xlsx"

# === READ TPM DATA ===
cat("[INFO] Reading TPM file:", tpm_file, "\n")
tpm_data <- read.table(tpm_file, header = TRUE, sep = "\t", check.names = FALSE)
cat("[INFO] TPM file loaded with", nrow(tpm_data), "rows and", ncol(tpm_data), "columns.\n")

# Filter for exon 1
tpm_exon1 <- tpm_data %>%
  filter(grepl("-E1$", Geneid))

# Extract PfHB3 gene IDs (remove exon info)
tpm_exon1 <- tpm_exon1 %>%
  mutate(GeneID_clean = sub("exon_(PfHB3_\\d+)-.*", "\\1", Geneid))

# === READ VAR GENE LIST ===
cat("[INFO] Reading var gene list:", var_gene_file, "\n")
var_genes <- read.table(var_gene_file, header = FALSE)
colnames(var_genes) <- c("GeneID")
cat("[INFO] Var gene list loaded with", nrow(var_genes), "entries.\n")

var_gene_ids <- var_genes$GeneID

# Filter TPM table for matching var genes
var_exon1_tpm <- tpm_exon1 %>%
  filter(GeneID_clean %in% var_gene_ids)

# Select and rename relevant columns
final_table <- var_exon1_tpm %>%
  select(GeneID_clean, `pipeline_output.bam.combined.bam`) %>%
  rename(GeneID = GeneID_clean, TPM = `pipeline_output.bam.combined.bam`)

# Write to Excel
write_xlsx(final_table, output_file)

cat("[INFO] Exported", nrow(final_table), "var genes (exon 1) to", output_file, "\n")


