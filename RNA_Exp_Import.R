# This code is used to import the RNA expression from the TCGA UCEC Data set from the GDC (https://gdc.cancer.gov/)

library("TCGAbiolinks")
library("limma")
library("edgeR")
library("dplyr")
library('SummarizedExperiment')

### Should only need to run once. 
query_TCGA = GDCquery(
  project = "TCGA-UCEC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts", 
  sample.type=c("Primary Tumor"),
  data.type="Gene Expression Quantification")

ucec_res = getResults(query_TCGA) # make results as table
GDCdownload(query = query_TCGA)
tcga_data_RNA = GDCprepare(query_TCGA)

saveRDS(object = tcga_data_RNA,
        file = "ucec_starCounts.RDS")
