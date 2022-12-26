library(pheatmap)
#tpm 
#delet the gene which expression is too low
breast_gtex_logtpm <- log2(breast_gtex_tpm+1)
breast_gtex_logtpm <- as.matrix(breast_gtex_logtpm)
gene_notzero_num_tpm <- apply(breast_gtex_logtpm[,image_exp_clu_median], 1, FUN = function(x){length(which(x!=0))})
gene_notoolow_num_tpm <- apply(breast_gtex_logtpm[,image_exp_clu_median],1,FUN = function(x){length(which(x>1e-08))})
result_feature_kmeans_median_tpm <- read.csv('gene_analysis.csv')
#tidy the result
result_feature_kmeans_median_tpm <- t(na.omit(t(result_feature_kmeans_median_tpm)))
result_feature_kmeans_median_tpm <- t(result_feature_kmeans_median_tpm)
result_feature_kmeans_median_tpm <- as.data.frame(result_feature_kmeans_median_tpm)
#tpval，chpval
chpval <- 1-pchisq(result_feature_kmeans_median_tpm[,3],4-2)
tpval <- 2*(1-pt(abs(result_feature_kmeans_median_tpm[,2]),length(group_km_median)-4-1))
result_feature_kmeans_median_tpm$tpval <- tpval
result_feature_kmeans_median_tpm$chpval <- chpval
colnames(result_feature_kmeans_median_tpm)[c(1,2,3)] <- c("cluster","tsatstic","chstastic")
logfold_change <- c()
for (i in c(1:dim(result_feature_kmeans_median_tpm)[1])) {
  
  logfold_change[i] <- caculate_fold_change(breast_gtex_logtpm,
                                            ind1 = match_names[match(km_median_df$sample[which(as.numeric(km_median_df$cluster)==result_feature_kmeans_median_tpm$cluster[i])],match_names$subject_id),]$geneExp.sample.ID,
                                            ind2 = setdiff(image_exp_clu_median,match_names[match(km_median_df$sample[which(as.numeric(km_median_df$cluster)==result_feature_kmeans_median_tpm$cluster[i])],match_names$subject_id),]$geneExp.sample.ID),
                                            index = rownames(result_feature_kmeans_median_tpm)[i])
  
}
result_feature_kmeans_median_tpm$logfoldchange <- logfold_change
result_feature_kmeans_median_tpm <- result_feature_kmeans_median_tpm[order(result_feature_kmeans_median_tpm$logfoldchange,decreasing = T),]
median_value <- unlist(apply(breast_gtex_logtpm[rownames(result_feature_kmeans_median_tpm),image_exp_clu_median],1,median))
result_feature_kmeans_median_tpm$base <- median_value

#select the gene which tpval less than 0.01，chipval better than 0.1
obvious_feature_result_median_tpm <- subset(result_feature_kmeans_median_tpm,chpval>0.1&tpval<0.01)
obvious_feature_result_median_tpm <- arrange(obvious_feature_result_median_tpm,-obvious_feature_result_median_tpm$tsatstic)
obvious_feature_result_median_tpm <- obvious_feature_result_median_tpm[order(obvious_feature_result_median_tpm$cluster),]
obvious_feature_result_median_tpm$chpval <- p.adjust(obvious_feature_result_median_tpm$chpval, method = "BH")
obvious_feature_result_median_tpm$tpval <- p.adjust(obvious_feature_result_median_tpm$tpval, method = "BH")
obvious_lofC_tpm <- subset(obvious_feature_result_median_tpm,logfoldchange>=0.5)


#colormap
gene_km_median_list_tpm <- c(rownames(subset(obvious_feature_result_median_tpm,cluster==1))[c(1:25)],
                             rownames(subset(obvious_feature_result_median_tpm,cluster==2))[c(1:25)],
                             rownames(subset(obvious_feature_result_median_tpm,cluster==3))[c(1:25)],
                             rownames(subset(obvious_feature_result_median_tpm,cluster==4))[c(1:25)]
)

gene_annotate_tpm <- data.frame(cluster=factor(rep(c("1","2","3","4"),each=25)),row.names = gene_km_median_list_tpm)

rownames(km_median_df) <- match_names$geneExp.sample.ID[match(km_median_df$sample,match_names$subject_id)]
color_matrix_median_tpm <- breast_gtex_logtpm[gene_km_median_list_tpm,rownames(km_median_df)]
for (i in c(1:dim(color_matrix_median_tpm)[1])) {
  
  c = median(color_matrix_median_tpm[i,])
  for (j in c(1:dim(color_matrix_median_tpm)[2])) {
    
    color_matrix_median_tpm[i,j] <- ifelse(color_matrix_median_tpm[i,j]>c,1,0)
    
  }
  
}

map_color <- list(gene_belong=c(cluster_1="#E41A1C",cluster_2="#FFFF33",cluster_3="#4DAF4A",cluster_4="#FF7F00"),
                  cluster=c("1"="#8DD3C7","2"="#FDB462","3"="#FFFFB3","4"="#FB8072"))

pheatmap(color_matrix_median_tpm,
         cluster_cols = F,cluster_rows = F,show_colnames = F,scale = "none",
         color = c("green","red"),
         annotation_col = select(km_median_df,"cluster"),
         annotation_row = gene_annotate_tpm,
         legend_breaks = c(1,0),
         legend_labels = c("high","low"),
         fontsize_row = 5,
         show_rownames = F)


box_tpm_df <- data.frame(gene=breast_gtex_logtpm["MT-TM",image_exp_clu_median],samples = image_exp_clu_median,cluster=km_median_df$cluster)
drawbox <- function(genenames,clustern){
  
  box_tpm_df$gene <- breast_gtex_logtpm[genenames,image_exp_clu_median]
  boxplot(gene~cluster,data = box_tpm_df,col=brewer.pal(4,"Set3"),main=paste0("Obvious in cluster",as.character(clustern)),
          ylab = genenames)
  
  
  
}
par(mfrow=c(3,3))
for (i in c(1:9)) {
  
  
  drawbox(rownames(subset(obvious_feature_result_median_tpm,cluster==2))[i],2)
  
}

