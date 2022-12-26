library(tidyr)
library(dplyr)

#select_feature
image_feature_select <- read.csv("features.csv")
image_feature_select <- column_to_rownames(image_feature_select,var = "subject_id")

for (i in grep(colnames(image_feature_select),pattern = "_hist",fixed = T)){
  image_feature_select[i] <- unlist(lapply(image_feature_select[i],FUN = function(x){x <- gsub("[","",x,fixed = T)}))
  image_feature_select[i] <- unlist(lapply(image_feature_select[i],FUN = function(x){x <- gsub("]","",x,fixed = T)}))
}

for (i in colnames(image_feature_select)[grep(colnames(image_feature_select),pattern = "_hist",fixed = T)]){
  a <- lapply(image_feature_select[,i],FUN = function(x){unlist(strsplit(x,",",fixed = T))})
  a <- lapply(a,as.numeric)
  a <- lapply(a,FUN = function(x){x/sum(x)})
  a <- lapply(a,as.character)
  a <- lapply(a, function(x){a <- paste(x,collapse = ",")})
  image_feature_select[,i] <- unlist(a)
  image_feature_select <- separate(image_feature_select,i,into=paste0(i,c(1:10)),sep = ",") 
}
for (i in colnames(image_feature_select)[grep(colnames(image_feature_select),pattern = "_hist",fixed = T)]){
  image_feature_select[,i] <- as.numeric(image_feature_select[,i]) 
}




#select median,mean,stddev，delet all_zero and duplicated_feature
feature_notzero_num <- apply(image_feature_select, 2, FUN = function(x){length(which(x!=0))})
image_feature_select_median <- image_feature_select[,which(feature_notzero_num>10)]
image_feature_select_median <- select(image_feature_select_median,ends_with(c("_median","_mean","_stddev")))
image_feature_select_median <- select(image_feature_select_median,-contains("Flatness"))
image_feature_select_median <- select(image_feature_select_median,-contains(c("NumberOfPixelsOnBorder")))
image_feature_select_median <- select(image_feature_select_median,-contains(c("voronoi","delaunay","mst")))
image_feature_select_median <- select(image_feature_select_median,-contains(c("density_distance_for_neighbors")))


#Feature category
feature_message <- fromJSON(file = "features.json")
feature_message <- sapply(feature_message,function(x){
  feature_name <- names(x)
  feature_type <- x
})
feature_categories <- data.frame(feature_name=names(feature_message),feature_type=feature_message)
b <- grep("_hist",feature_categories$feature_name,fixed = T)
c <- data.frame(feature_name=paste0(rep(feature_categories$feature_name[b],each=10),rep(c("1","2","3","4","5","6","7","8","9","10"),times=10)),
                feature_type=rep(feature_categories$feature_type[b],each=10))
feature_categories <- feature_categories[-b,]%>%
  rbind(c)
rownames(feature_categories) <- feature_categories$feature_name
feature_categories_median <- feature_categories[match(colnames(image_feature_select_median),feature_categories$feature_name),]
feature_categories_median <- select(feature_categories_median,"feature_type")
feature_categories_median <- arrange(feature_categories_median,feature_categories_median[,1])
