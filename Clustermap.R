library(ggplot2)
library(RColorBrewer)
library(pheatmap)
#comfirm the number of clusters
k_sse3 <- c()
k_sse3_round <- c()
for (i in c(2:15)) {
  
  for(j in c(1:5)){
    
    k_sse3_round[j] <- kmeans(scale(image_feature_select_median),centers = i)[[5]]
  }
  
  k_sse3[i-1] <- mean(k_sse3_round)
  
}
plot(k_sse3)
lines(k_sse3)
plot(-diff(x),type='b')

#UMAP
library(umap)
k_list <- list()
for(i in c(2:6)){
  k_list[[i]] <- kmeans(customer, centers = i)$cluster 
}
image.umap <- umap::umap(scale(image_feature_select_median))
umap_data <- as.data.frame(image.umap$layout, row.names = rownames(customer))
colnames(umap_data) <- c("UMAP1","UMAP2")
umap_data$group <- as.character(k_list[[i]])#select the number to show
ggplot(data = umap_data, mapping = aes(UMAP1,UMAP2,color=group))+
  geom_point(size = 2.0, shape = 16)+
  scale_fill_manual(values=c("red","blue","green","yellow","purple","gray"))+
  ggtitle("UMAP")


##MLP_K-Means:new_result
km_median_df <- read.csv('./result.csv')
km_median_df <- select(km_median_df,"cluster")
km_median_df$cluster <- km_median_df$cluster+1
km_median_df$sample <- rownames(image_feature_select_median)
rownames(km_median_df) <- km_median_df$sample


km_median_df <- km_median_df[order(km_median_df$cluster),]
km_median_df$cluster <- as.factor(km_median_df$cluster)
ann_color <- list(feature_type=c(distance_feature="#E41A1C",intensity_feature="#FFFF33",shape_feature="#4DAF4A",size_feature="#FF7F00"),
                  cluster=c("1"="#8DD3C7","2"="#FDB462","3"="#FFFFB3","4"="#FB8072"))
pheatmap(image_feature_select_median[km_median_df$sample,rownames(feature_categories_median)],
         cluster_cols = F,cluster_rows = F,show_rownames = F,scale = 'column',
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(101),
         annotation_row = select(km_median_df,"cluster"),annotation_col = feature_categories_median,
         annotation_colors = ann_color,
         breaks=seq(-2,2,by=0.04))

image_exp_clu_median <- match_names$geneExp.sample.ID[match(rownames(km_median_df),match_names$subject_id)]
group_km_median <- km_median_df$cluster
