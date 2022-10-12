
# This of code was used to extract the candidate gene of interest expression data in cell line AN3-CA from the RNA cell line consensus dataset downloaded from cBioportal (https://www.cbioportal.org/).

# Attach required library

library(dplyr)

# Import downloaded RNA cell line dataset. 

Cell_line_dataset=read.delim("rna_celline.tsv")

# Subset RNA cell line dataset to include only gene expression data from cell line AN3-CA.

AN <- subset(Cell_line_dataset, Cell.line == "AN3-CA", select = c("Gene.name", "nTPM"))

# Subset AN3-CA extracted dataset to include only candidate genes of interest. 

ANExpress <- subset(AN, Gene.name %in% c("SLCO1B3", "SALL3", "SLC6A3", "LPCAT1", "NXPH1", "STK11", "SPINT4", "PTGIS", "PMS2", "MAPK3", "SPN", "QPRT", "C16orf54", "MAZ", "PRRT2", "MVP", "CDIPT",
                                         "SEZ6L2", "ASPHD1", "KCTD13", "TMEM219", "TAOK2", "HIRIP3", "DOC2A", "TLCD3B", "ALDOA", "PPP4C", "TBX6", "YPEL3", "GDPD3", "MSH2", "ZNF185", "WFDC3", "ASXL3", "CBARP", "CYLC2",
                                         "MSH6", "FBXO11", "TAF4", "ITPR1", "BRCA1", "GLDC", "TERT", "NAT1", "ZNF320", "ATP2B2", "CDKN1C", "C7orf50", "ERCC2", "DAB2IP", "RANBP9", "SMARCA2", "SKI", "JRKL", "ERBB4", "CD36",
                                         "ZSCAN5A", "UPK3B", "ARHGEF10", "AUTS2", "TBX1", "VCX3A", "DIP2B", "NEIL2", "LAMA5", "SOD2", "MAGI1", "LARGE1", "SLCO1B1", "BARD1", "MUTYH", "AKT1", "XRCC1"))


