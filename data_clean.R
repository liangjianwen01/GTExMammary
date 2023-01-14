#Expression Matrix
# Readin GTEx gene read count after TPM
GTEx=read.table('TPM matrix path',header = TRUE,sep = '\t',skip = 2,check.names = F)
## Annotations
# SMTS Tissue Type, area from which the tissue sample was taken. 
# SMTSD Tissue Type, more specific detail of tissue type

b=read.table('annotation path',header = TRUE,sep = '\t',quote = '', check.names = F)
for (i in grep("(", b$SMTSD, fixed = TRUE)){
  
  b$SMTSD[i] <- strsplit(b$SMTSD[i], " (", fixed = TRUE)[[1]][1]
}

breast_gtex_tpm=GTEx[,colnames(GTEx) %in% b[b$SMTS=='Breast - Mammary Tissue',1]] 
rownames(breast_gtex_tpm)=GTEx[,1]

