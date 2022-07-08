args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]

source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")

## Format Expression file
expr = read.csv(file.path(input_dir, "EXPR.csv.gz")  , sep=";", stringsAsFactors=FALSE )
write.table( expr , file= file.path(output_dir, "EXPR.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=TRUE )

## Format clinical data
clin_original = read.csv( file.path(input_dir, "CLIN.txt"), stringsAsFactors=FALSE , sep="\t")

selected_cols <- c( "ID" , "Sex" , "AgeDGM" , "Histology" , "Tumor.stage" , "PFS" , "evtPFS" )
clin = cbind( clin_original[ , selected_cols] , "Lung" , "PD-1/PD-L1" , NA , NA, NA , NA , NA, NA , NA )
colnames(clin) = c( "patient" , "sex" , "age"  ,"histo"  , "stage" , "t.pfs" , "pfs" , "primary" , "drug_type" , "response.other.info" , "response"  , "t.os" , "os", "recist" , "dna" , "rna")

clin$rna = "tpm"
clin$response = Get_Response( data=clin )

clin$stage = ifelse( clin$stage %in% 3 , "III" , 
				ifelse( clin$stage %in% 4 , "IV" , 
				ifelse( clin$stage %in% c(0,1,2) , "I/II" , NA ) )) 
				

clin = clin[ , c("patient" , "sex" , "age" , "primary" , "histo" , "stage" , "response.other.info" , "recist" , "response" , "drug_type" , "dna" , "rna" , "t.pfs" , "pfs" , "t.os" , "os" ) ]
clin$patient = paste( "P" , clin$patient , sep="" )

clin_original$ID = paste( "P" , clin_original$ID , sep="" )
rownames(clin) <- clin$patient
clin = format_clin_data(clin_original, "ID", selected_cols, clin)

write.table( clin , file=file.path(output_dir, "CLIN.csv") , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )

