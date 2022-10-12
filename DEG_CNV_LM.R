## The aim of this code is to apply a function to perform filtration, normalization and perform a linear regression for each CNV comparing samples expression with or with a CNV deletion. 

# Create DGE list using matched RNA data in file "RNA_CNV_Match_nodup.R"

dge=DGEList(counts=assay(RNA),
            samples=colData(RNA),
            genes=rowData(RNA)$gene_name)

# Filtering the DGE object

keep <- filterByExpr(dge)
dge = dge[keep,,keep.lib.sizes=FALSE] 


## Normalize (Scaling and voom)
dge <- calcNormFactors(dge,method="TMM")
dge$samples$norm.factors
lcpm.normalised <- cpm(dge, log=TRUE)

v = voom(dge$counts,design,plot=TRUE)


# Function to perform the whole process

dge_list <- apply(CN[3:ncol(CN)], 2, function(i){
  group = factor(i)
  design <- model.matrix(~0+group) 
  colnames(design) <- c("del", "notDel")
  
  
  contr.matrix <- makeContrasts(
    delvsnot = del - notDel,
    levels=colnames(design)
  )
  
  fit <- lmFit(voom(dge$counts,design,plot=TRUE), design)
  fit <- contrasts.fit(fit, contrasts=contr.matrix)
  efit <- eBayes(fit)  
})


## Extraction of top 100 DEGs

dge_list_topTable <- lapply(dge_list, function(i){
  topTable(i, n=100, genelist=dge$genes)
} )

dge_summary <- lapply(dge_list, function(i){
  summary(decideTests(i))
})
dge_list_topTable[[]]
