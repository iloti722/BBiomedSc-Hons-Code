# The code in this file is used to match the patient samples in the CNV ateration dataset with patient samples in the RNA dataset, removing duplicates. 

#Import files downloaded from cbioportal and R file "RNA_Exp_Import"

RNA <- readRDS("ucec_starCounts.RDS")
CN <- read_xlsx("TCGA_CNV_newdownload.xlsx")

# Sort CN Data
# Altering CNV data to only be 0,1 (1 = Diploid, Gain, Amp, 0 = Shallow Del, Deep Del)

CN <- CN %>% mutate(across(3:70, ~ifelse(.x >=0, 1, 0)))


## Subset RNA
tmp<-RNA[,colData(RNA)$patient %in% gsub("-\\w+$", "", CN$SAMPLE_ID)]

x<-as.data.frame(colData(tmp))%>% group_by(patient) %>% summarise(n=n()) %>% filter(n>1) %>% pull(patient)
duplicated_samples =as.data.frame(colData(tmp)) %>% filter(patient %in% x) %>% arrange(patient)

# Removing duplicate samples in RNA
id2remove <- c('TCGA-BK-A26L-01C-04R-A277-07', 'TCGA-BK-A26L-01A-11R-A16F-07', 'TCGA-BK-A0CA-01B-02R-A277-07', 'TCGA-BK-A0CA-01A-21R-A118-07', 'TCGA-BK-A0CC-01B-04R-A277-07', 'TCGA-BK-A0CC-01A-21R-A16W-07', 'TCGA-BK-A139-01A-11R-A118-07', 'TCGA-BK-A139-01C-08R-A277-07') ## List of ids to remove
RNA<-tmp[,!rownames(colData(tmp)) %in%id2remove] ## remove sample matching id2remove


## match CNV row (Sample_ID) to RNA patient IDs 
CN<-CN[ match(colData(RNA)$patient, gsub("-\\w+$", "", CN$SAMPLE_ID)),]
