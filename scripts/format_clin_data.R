# wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190266/matrix/GSE190266_series_matrix.txt.gz' \
# -O ~/Code/bhklab/ICBCuration/download/ICB_Fumet1/GSE190266_series_matrix.txt.gz

library(data.table)
library(stringr)
library(GEOquery)
library(tibble)

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

gunzip(file.path(input_dir, "GSE190266_series_matrix.txt.gz"))
clin <- getGEO(filename = file.path(input_dir, "GSE190266_series_matrix.txt"), destdir = input_dir)
clin <- pData(clin)
colnames(clin)[colnames(clin) == "title"] <- "patient"
clin$patient <- str_replace_all(clin$patient, "\\W", "_")
rownames(clin) <- clin$patientid

# TODO: Format the clinical data with common columns etc.
selected_cols <- c("patient", "tissue:ch1", "pfs_time (6 months):ch1", "pfs_evt (6 months):ch1")
remaining_cols <- colnames(clin)[!colnames(clin) %in% selected_cols]
clin <- clin[, c(selected_cols, remaining_cols)]
colnames(clin)[colnames(clin) %in% selected_cols] <- c("patient", "tissueid", "t.pfs", "pfs")
clin$tissueid <- "Lung"
clin <- add_column(clin, response = NA, response.other.info = NA, recist = NA, .after = "tissueid")
clin$t.pfs <- as.numeric(clin$t.pfs)
clin$pfs <- as.numeric(clin$pfs)
clin$response <- Get_Response(clin)

clin <- add_column(
  clin,
  sex = NA,
  age = NA,
  primary = "Lung",
  stage = NA,
  histo = NA,
  treatmentid = "",
  drug_type = "anti-PD-1",
  dna = NA,
  rna = "tpm",
  t.os = NA,
  os = NA,
  survival_unit = "month",
  survival_type = NA,
  .after = "patient"
)

clin <- clin[, c(
  "patient", "sex", "age", "primary", "treatmentid", "tissueid", "histo", "stage", "response.other.info", "recist", 
  "response", "drug_type", "dna", "rna", "t.pfs", "pfs", "t.os", "os", "survival_unit", "survival_type", 
  remaining_cols
)]
# colnames(clin)[colnames(clin) %in% c("t.pfs", "pfs", "t.os", "os")] <- c("survival_time_pfs", "event_occurred_pfs", "survival_time_os", "event_occurred_os")
rownames(clin) <- clin$patient
write.table(clin, file = file.path(output_dir, "CLIN.csv"), quote = FALSE, sep = ";", col.names = TRUE, row.names = TRUE)

case <- as.data.frame(cbind(clin$patient, rep(0, length(clin$patient)), rep(0, length(clin$patient)), rep(1, length(clin$patient))))
colnames(case) <- c("patient", "snv", "cna", "expr")
rownames(case) <- clin$patient
write.table(case, file = file.path(output_dir, "cased_sequenced.csv"), quote = FALSE, sep = ";", col.names = TRUE, row.names = TRUE)
