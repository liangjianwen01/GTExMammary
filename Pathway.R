library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

diff1 <- subset(obvious_feature_result_median_tpm,cluster==1)
diff2 <- subset(obvious_feature_result_median_tpm,cluster==2)
diff3 <- subset(obvious_feature_result_median_tpm,cluster==3)
diff4 <- subset(obvious_feature_result_median_tpm,cluster==4)


gene.df1 <- bitr(geneID=rownames(diff1) ,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene.df2 <- bitr(geneID=rownames(diff2) ,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
gene.df3 <- bitr(geneID=rownames(diff3) ,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)                
gene.df4 <- bitr(geneID=rownames(diff4) ,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)

gene1 <- gene.df1$ENTREZID
gene2 <- gene.df2$ENTREZID
gene3 <- gene.df3$ENTREZID
gene4 <- gene.df4$ENTREZID

ego_ALL1 <- enrichGO(gene = gene1,
                     OrgDb=org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
ego_ALL2 <- enrichGO(gene = gene2,
                     OrgDb=org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
ego_ALL3 <- enrichGO(gene = gene3,
                     OrgDb=org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
ego_ALL4 <- enrichGO(gene = gene4,
                     OrgDb=org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
ego_ALL1 <- as.data.frame(ego_ALL1)
ego_ALL2 <- as.data.frame(ego_ALL2)
ego_ALL3 <- as.data.frame(ego_ALL3)
ego_ALL4 <- as.data.frame(ego_ALL4)

ego_ALL2 <- ego_ALL2[order(ego_ALL2$p.adjust),]
go_enrich_df2 <- data.frame(
  ID=ego_ALL2$ID[c(1:10)],
  Description=ego_ALL2$Description[c(1:10)],
  GeneNumber=ego_ALL2$Count[c(1:10)],
  padj <- ego_ALL2$p.adjust[c(1:10)],
  type=factor(ego_ALL2$ONTOLOGY[c(1:10)], 
              labels=c("biological process", "cellular component","molecular function" )))

ego_ALL3 <- ego_ALL3[order(ego_ALL3$p.adjust),]
go_enrich_df3 <- data.frame(
  ID=ego_ALL3$ID[c(1:10)],
  Description=ego_ALL3$Description[c(1:10)],
  GeneNumber=ego_ALL3$Count[c(1:10)],
  padj <- ego_ALL3$p.adjust[c(1:10)],
  type=factor(ego_ALL3$ONTOLOGY[c(1:10)], 
              labels=c("biological process", "cellular component","molecular function" )))

ego_ALL4 <- ego_ALL4[order(ego_ALL4$p.adjust),]
go_enrich_df4 <- data.frame(
  ID=ego_ALL4$ID[c(1:8)],
  Description=ego_ALL4$Description[c(1:8)],
  GeneNumber=ego_ALL4$Count[c(1:8)],
  padj <- ego_ALL4$p.adjust[c(1:8)],
  type=factor(ego_ALL4$ONTOLOGY[c(1:8)], 
              labels=c("biological process", "molecular function" )))

for(i in 1:nrow(go_enrich_df2)){
  description_splite=strsplit(go_enrich_df2$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:6],collapse = " ") 
  go_enrich_df2$Description[i]=description_collapse
  go_enrich_df2$Description=gsub(pattern = "NA","",go_enrich_df2$Description)
}
for(i in 1:nrow(go_enrich_df3)){
  description_splite=strsplit(go_enrich_df3$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:6],collapse = " ") 
  go_enrich_df3$Description[i]=description_collapse
  go_enrich_df3$Description=gsub(pattern = "NA","",go_enrich_df3$Description)
}
for(i in 1:nrow(go_enrich_df4)){
  description_splite=strsplit(go_enrich_df4$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:4],collapse = " ") 
  go_enrich_df4$Description[i]=description_collapse
  go_enrich_df4$Description=gsub(pattern = "NA","",go_enrich_df4$Description)
}


go_enrich_df2 <- go_enrich_df2[order(go_enrich_df2$type),]
rownames(go_enrich_df2) <- as.character(c(1:10))
go_enrich_df2$type_order=factor(rev(as.integer(rownames(go_enrich_df2))),labels=rev(go_enrich_df2$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
go_enrich_df2[,4] <- -log10(go_enrich_df2$padj....ego_ALL2.p.adjust.c.1.10..)
colnames(go_enrich_df2)[4] <- 'Nlog10_Pvalue'
ggplot(data=go_enrich_df2, aes(x=type_order,y=Nlog10_Pvalue, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("-log10(P-value)") + 
  labs(title = "The Most Enriched GO Terms in Cluster2")+
  theme(axis.text = element_text(size = 15))+
  theme(plot.title = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  theme(legend.position="none")

go_enrich_df3 <- go_enrich_df3[order(go_enrich_df3$type),]
rownames(go_enrich_df3) <- as.character(c(1:10))
go_enrich_df3$type_order=factor(rev(as.integer(rownames(go_enrich_df3))),labels=rev(go_enrich_df3$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
go_enrich_df3[,4] <- -log10(go_enrich_df3$padj....ego_ALL3.p.adjust.c.1.10..)
colnames(go_enrich_df3)[4] <- 'Nlog10_Pvalue'
ggplot(data=go_enrich_df3, aes(x=type_order,y=Nlog10_Pvalue, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("-log10(P-value)") + 
  labs(title = "The Most Enriched GO Terms in Cluster3")+
  theme(axis.text = element_text(size = 12))+
  theme(plot.title = element_text(size = 15))+
  theme(axis.title = element_text(size = 15))+
  theme(legend.position="none")

go_enrich_df4 <- go_enrich_df4[order(go_enrich_df4$type),]
rownames(go_enrich_df4) <- as.character(c(1:8))
go_enrich_df4$type_order=factor(rev(as.integer(rownames(go_enrich_df4))),labels=rev(go_enrich_df4$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
go_enrich_df4[,4] <- -log10(go_enrich_df4$padj....ego_ALL4.p.adjust.c.1.8..)
colnames(go_enrich_df4)[4] <- 'Nlog10_Pvalue'
ggplot(data=go_enrich_df4, aes(x=type_order,y=Nlog10_Pvalue, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) +
  coord_flip() + 
  xlab("GO term") + 
  ylab("-log10(P-value)") + 
  labs(title = "The Most Enriched GO Terms in Cluster4")+
  theme(axis.text = element_text(size = 15))+
  theme(plot.title = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  theme(legend.position="none")


kk1 <- enrichKEGG(gene = gene1,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)
kk2 <- enrichKEGG(gene = gene2,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)
kk3 <- enrichKEGG(gene = gene3,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)
kk4 <- enrichKEGG(gene = gene4,keyType = "kegg",organism= "human", qvalueCutoff = 0.05, pvalueCutoff=0.05)

hh1 <- as.data.frame(kk1)

hh2 <- as.data.frame(kk2)

hh3 <- as.data.frame(kk3)

hh4 <- as.data.frame(kk4)


rownames(hh1) <- 1:nrow(hh1)
hh1$p.adjust <- -log10(hh1$p.adjust)
hh1$order=factor(rev(as.integer(rownames(hh1))),labels = rev(hh1$Description))
ggplot(hh1,aes(y=order,x=p.adjust,fill='red'))+
  geom_bar(stat = "identity",width=0.7)+
  labs(title = "KEGG Pathways Enrichment in Cluster1",
       x = "-log10(P-value)", 
       y = "Pathways")+
  theme(axis.text = element_text(size = 15))+
  theme(plot.title = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  theme(legend.position="none")

rownames(hh2) <- 1:nrow(hh2)
hh2$p.adjust <- -log10(hh2$p.adjust)
hh2$order=factor(rev(as.integer(rownames(hh2))),labels = rev(hh2$Description))
ggplot(hh2,aes(y=order,x=p.adjust,fill='red'))+
  geom_bar(stat = "identity",width=0.7)+
  labs(title = "KEGG Pathways Enrichment in Cluster2",
       x = "-log10(P-value)", 
       y = "Pathways")+
  theme(axis.text = element_text(size = 15))+
  theme(plot.title = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  theme(legend.position="none")

rownames(hh3) <- 1:nrow(hh3)
hh3$p.adjust <- -log10(hh3$p.adjust)
hh3$order=factor(rev(as.integer(rownames(hh3))),labels = rev(hh3$Description))
ggplot(hh3,aes(y=order,x=p.adjust,fill='red'))+
  geom_bar(stat = "identity",width=0.7)+
  labs(title = "KEGG Pathways Enrichment in Cluster3",
       x = "-log10(P-value)", 
       y = "Pathways")+
  theme(axis.text = element_text(size = 15))+
  theme(plot.title = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  theme(legend.position="none")

rownames(hh4) <- 1:nrow(hh4)
hh4$p.adjust <- -log10(hh4$p.adjust)
hh4$order=factor(rev(as.integer(rownames(hh4))),labels = rev(hh4$Description))
ggplot(hh4,aes(y=order,x=p.adjust,fill='red'))+
  geom_bar(stat = "identity",width=0.7)+
  labs(title = "KEGG Pathways Enrichment in Cluster4",
       x = "-log10(P-value)", 
       y = "Pathways")+
  theme(axis.text = element_text(size = 15))+
  theme(plot.title = element_text(size = 18))+
  theme(axis.title = element_text(size = 18))+
  theme(legend.position="none")
