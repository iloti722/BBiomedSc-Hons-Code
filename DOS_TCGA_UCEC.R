# The aim of this code is to calculate whether the RNA expression profiles of different copy number alterations are significantly different, indicating dosage sensitivity.

# Install packages

library(cBioPortalData)
library(dplyr)
library(ggplot2)

# Set up URL for extraction

cbio <- cBioPortal(
  hostname = "www.cbioportal.org",
  protocol = "https",
  api. = "/api/api-docs")

# Select correct study (TCGA UCEC)

cbio <- cBioPortal()
studies <- getStudies(cbio)
studies %>% filter(grepl("Endo", name)) %>% select(name, allSampleCount,studyId) %>% as.data.frame()
studyId='ucec_tcga_pan_can_atlas_2018'      

## list the molecular data available

mols <- molecularProfiles(cbio,studyId)
mols[["molecularProfileId"]]

## Extract CNV alteration data from TCGA UCEC from cbioportal
x<-getDataByGenes(
  api = cbio,
  studyId = studyId,
  genes = c("SLCO1B3", "SALL3", "SLC6A3", "LPCAT1", "NXPH1", "STK11", "SPINT4", "PTGIS", "PMS2", "MAPK3", "SPN", "QPRT", "C16orf54", "MAZ", "PRRT2", "MVP", "CDIPT",
            "SEZ6L2", "ASPHD1", "KCTD13", "TMEM219", "TAOK2", "HIRIP3", "DOC2A", "TLCD3B", "ALDOA", "PPP4C", "TBX6", "YPEL3", "GDPD3", "MSH2", "ZNF185", "WFDC3", "ASXL3", "CBARP", "CYLC2",
            "MSH6", "FBXO11", "TAF4", "ITPR1", "BRCA1", "GLDC", "TERT", "NAT1", "ZNF320", "ATP2B2", "CDKN1C", "C7orf50", "ERCC2", "DAB2IP", "RANBP9", "SMARCA2", "SKI", "JRKL", "ERBB4", "CD36",
            "ZSCAN5A", "UPK3B", "ARHGEF10", "AUTS2", "TBX1", "VCX3A", "DIP2B", "NEIL2", "LAMA5", "SOD2", "MAGI1", "LARGE1", "SLCO1B1", "BARD1", "MUTYH", "AKT1", "XRCC1"),
  by = "hugoGeneSymbol",
  molecularProfileIds = "ucec_tcga_pan_can_atlas_2018_gistic"
)

## Extract RNA expression data from the TCGA UCEC cbioportal
y<-getDataByGenes(
  api =  cbio,
  studyId = studyId,
  genes = c("SLCO1B3", "SALL3", "SLC6A3", "LPCAT1", "NXPH1", "STK11", "SPINT4", "PTGIS", "PMS2", "MAPK3", "SPN", "QPRT", "C16orf54", "MAZ", "PRRT2", "MVP", "CDIPT",
            "SEZ6L2", "ASPHD1", "KCTD13", "TMEM219", "TAOK2", "HIRIP3", "DOC2A", "TLCD3B", "ALDOA", "PPP4C", "TBX6", "YPEL3", "GDPD3", "MSH2", "ZNF185", "WFDC3", "ASXL3", "CBARP", "CYLC2",
            "MSH6", "FBXO11", "TAF4", "ITPR1", "BRCA1", "GLDC", "TERT", "NAT1", "ZNF320", "ATP2B2", "CDKN1C", "C7orf50", "ERCC2", "DAB2IP", "RANBP9", "SMARCA2", "SKI", "JRKL", "ERBB4", "CD36",
            "ZSCAN5A", "UPK3B", "ARHGEF10", "AUTS2", "TBX1", "VCX3A", "DIP2B", "NEIL2", "LAMA5", "SOD2", "MAGI1", "LARGE1", "SLCO1B1", "BARD1", "MUTYH", "AKT1", "XRCC1"),
  by = "hugoGeneSymbol",
  molecularProfileIds = "ucec_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
)


## create patient/sample and gene identifier.
x1 <- x$ucec_tcga_pan_can_atlas_2018_gistic %>% mutate(id = paste0(sampleId,hugoGeneSymbol))
y1 <- y$ucec_tcga_pan_can_atlas_2018_rna_seq_v2_mrna %>% mutate(id = paste0(sampleId,hugoGeneSymbol))


## join data, filter, plot
left_join(x1,y1, by=c("id"="id"))%>%
  filter(!is.na(value.x), !is.na(value.y))%>%ggplot(aes(x=value.x, y=value.y, group=value.x)) +geom_boxplot() +
  geom_point()+facet_wrap(hugoGeneSymbol.x~., scales = 'free_x')

#~~~~~~~~~~~~~~~~~~~~NEW CODE~~~~~~~~~~~~~~~~~~~~#

## join data/filter out missing (x and y must be same length for cor.test) and test for association between copy number alteration and RNA expresssion.

left_join(x1,y1, by=c("id"="id"))%>%filter(!is.na(value.x), !is.na(value.y)) %>%
  group_by(hugoGeneSymbol.x) %>% summarise(p=cor.test(value.x,value.y)$p.value)

?cor.test
