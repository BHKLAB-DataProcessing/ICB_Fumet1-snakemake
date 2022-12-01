from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]
filename = config["filename"]

rule get_multiassay:
    input:
        S3.remote(prefix + "scripts/get_multiassay.R"),
        prefix + "processed/EXPR_gene_tpm.rds",
        prefix + "processed/EXPR_gene_counts.rds",
        prefix + "processed/EXPR_isoform_tpm.rds",
        prefix + "processed/EXPR_isoform_counts.rds"
    output:
        prefix + filename
    resources:
        mem_mb = 3000,
        disk_mb = 3000
    shell:
        """
        Rscript {prefix}scripts/get_multiassay.R \
        {prefix}processed \
        {prefix} \
        ICB_Limagne1 \
        EXPR_gene_tpm:EXPR_gene_counts:EXPR_isoform_tpm:EXPR_isoform_counts    
        """

rule get_exp_se:
    input:
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData"),
        S3.remote(prefix + "scripts/get_exp_se.R"),
        prefix + "processed/EXPR_gene_tpm.csv",
        prefix + "processed/EXPR_gene_counts.csv",
        prefix + "processed/EXPR_isoform_tpm.csv",
        prefix + "processed/EXPR_isoform_counts.csv",
        prefix + "processed/CLIN.rds",
        prefix + "processed/cased_sequenced.rds"
    output:
        prefix + "processed/EXPR_gene_tpm.rds",
        prefix + "processed/EXPR_gene_counts.rds",
        prefix + "processed/EXPR_isoform_tpm.rds",
        prefix + "processed/EXPR_isoform_counts.rds"
    shell:
        """
        Rscript {prefix}scripts/get_exp_se.R \
        {prefix}processed \
        Fumet1 \
        FALSE \
        FALSE \
        TRUE \
        {prefix}annotation/Gencode.v40.annotation.RData        
        """

rule format_clin_cased_sequenced:
    input:
        S3.remote(prefix + "scripts/format_clin_cased_sequenced.R"),
        prefix + "processed/CLIN.csv",
        prefix + "processed/cased_sequenced.csv"
    output:
        prefix + "processed/cased_sequenced.rds",
        prefix + "processed/CLIN.rds"
    shell:
        """
        Rscript {prefix}scripts/format_clin_cased_sequenced.R \
        {prefix}processed
        """

rule download_scripts:
    output:
        S3.remote(prefix + "scripts/get_exp_se.R"),
        S3.remote(prefix + "scripts/get_multiassay.R"),
        S3.remote(prefix + "scripts/format_clin_cased_sequenced.R")
    shell:
        """
        wget https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/code/get_exp_se.R -O {prefix}scripts/get_exp_se.R
        wget https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/code/get_multiassay.R -O {prefix}scripts/get_multiassay.R
        wget https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/code/format_clin_cased_sequenced.R -O {prefix}scripts/format_clin_cased_sequenced.R
        """

rule format_expr_data:
    input:
        S3.remote(prefix + "download/expr_list.rds")
    output:
        prefix + "processed/EXPR_gene_tpm.csv",
        prefix + "processed/EXPR_gene_counts.csv",
        prefix + "processed/EXPR_isoform_tpm.csv",
        prefix + "processed/EXPR_isoform_counts.csv"
    shell:
        """
        Rscript scripts/format_expr_data.R \
        {prefix}download \
        {prefix}processed
        """

rule format_clin_data:
    input:
        S3.remote(prefix + "annotation/curation_tissue.csv"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "download/GSE190266_series_matrix.txt.gz")
    output:
        prefix + "processed/CLIN.csv",
        prefix + "processed/cased_sequenced.csv"
    shell:
        """
        Rscript scripts/format_clin_data.R \
        {prefix}download \
        {prefix}processed \
        {prefix}annotation \
        """

rule compile_expr_data:
    input:
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData"),
        S3.remote(prefix + "download/ICB_Fumet1_RNAseq.zip"),
        S3.remote(prefix + "download/filereport_read_run_PRJNA786567_tsv.txt")
    output:
        S3.remote(prefix + "download/expr_list.rds"),
    shell:
        """
        Rscript scripts/compile_expr_data.R \
        {prefix}download \
        {prefix}annotation \
        """

rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData"),
        S3.remote(prefix + "annotation/curation_tissue.csv"),
        S3.remote(prefix + "annotation/curation_drug.csv")
    shell:
        """
        wget https://github.com/BHKLAB-DataProcessing/Annotations/blob/master/Gencode.v40.annotation.RData?raw=true -O {prefix}annotation/Gencode.v40.annotation.RData 
        wget https://github.com/BHKLAB-DataProcessing/ICB_Common/raw/main/data/curation_drug.csv -O {prefix}annotation/curation_drug.csv
        wget https://github.com/BHKLAB-DataProcessing/ICB_Common/raw/main/data/curation_tissue.csv -O {prefix}annotation/curation_tissue.csv 
        """

rule download_data:
    output:
        S3.remote(prefix + "download/GSE190266_series_matrix.txt.gz"),
        S3.remote(prefix + "download/ICB_Fumet1_RNAseq.zip"),
        S3.remote(prefix + "download/filereport_read_run_PRJNA786567_tsv.txt")
    shell:
        """
        wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190266/matrix/GSE190266_series_matrix.txt.gz' \
            -O {prefix}download/GSE190266_series_matrix.txt.gz

        wget 'https://zenodo.org/record/7158251/files/ICB_Fumet1_RNAseq.zip?download=1' \
            -O {prefix}download/ICB_Fumet1_RNAseq.zip

        wget 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA786567&result=read_run&fields=sample_title&format=tsv&download=true&limit=0' \
            -O {prefix}download/filereport_read_run_PRJNA786567_tsv.txt
        """
