# wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190266/matrix/GSE190266_series_matrix.txt.gz' \
# -O ~/Code/bhklab/ICBCuration/download/ICB_Fumet1/GSE190266_series_matrix.txt.gz

library(data.table)
library(stringr)
library(GEOquery)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

gunzip(file.path(input_dir, "GSE190266_series_matrix.txt.gz"))
clin <- getGEO(filename = file.path(input_dir, "GSE190266_series_matrix.txt"), destdir = input_dir)
clin <- pData(clin)
colnames(clin)[colnames(clin) == "title"] <- "patientid"
clin$patientid <- str_replace_all(clin$patientid, "\\W", '_')
rownames(clin) <- clin$patientid

# TODO: Format the clinical data with common columns etc.


write.table(clin, file = file.path(output_dir, "ICB_Fumet1_metadata.tsv"), row.names = TRUE, col.names = TRUE, sep = "\t")
