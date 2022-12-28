# GTExMammary
Input the CSV table of image features to tidy.R for organizing your data and confirming the number of clusters.
Input the organized data to ./DeepClustering/main.py for clustering.
Input the clustering result to Clustermap.R for generating the heatmap.
Input the gene analysis result to gene_analysis.R to visualize the results.
Run Pathway.R to get KEGG and GO analysis results

# Summary
  Gene expression can be used to supplement information that cannot be obtained by visual inspection. Many studies have explored the relationship between phenotype and genotype in diseased tissues. 
  In this project, we comprehensively examined the relationship between nuclear features and gene expression in healthy breast cells. We used 456 healthy mammary gland samples from GTEx for analysis. We divided the 456 samples into 4 groups according to the features of the nuclei, and found highly expressed gene sets specific to each group. Finally, we performed pathway analysis for each gene set.
Related publication: Tian Mou. Jianwen Liang. et al. (2022) A comprehensive landscape of imaging feature associated RNA expression profiles in breast. Sensor.