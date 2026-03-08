# THIS SCRIPT WILL HELP YOU LEARN THE BASICS OF BULK RNA-SEQ DATA PROCESSING
# NOTE: The demo data used in this code is from GEO dataset GSE183947 (using only first 10 paired samples)

# Set your current working directory
setwd("C:/Users/user_name/folder")

# Set seed for reproducibility
set.seed(20)

#----------------------------------------------------------------------------------------------#
# PRINCIPAL COMPONENT ANALYSIS 
#----------------------------------------------------------------------------------------------#

# Load all required libraries
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(dplyr)
library(edgeR)

# Read the metadata:
metadata <- as.data.frame(read.csv('metadata.csv',header = T,row.names = 1))
dim(metadata) # this should be 20  3
head(metadata)
table(metadata$tissue)
# breast_tumor     normal_breast_tissue
#       10                   10
table(metadata$Donor)

# Read the raw counts:
counts <- as.data.frame(read.csv('raw_counts.csv',header = T,row.names = 1))
dim(counts) # 62700    20
counts[1:5,1:5]

# Keep the same samples in counts as in metadata:
length(intersect(colnames(counts),rownames(metadata))) # 20
setdiff(colnames(counts),rownames(metadata)) # 0

# Match the sample order
counts_sel <- counts[,rownames(metadata)]
dim(counts_sel) # 62700     20
all(colnames(counts_sel)==rownames(metadata))  # TRUE

# Filter out lowly expressed genes based on cpm:
cpm <- cpm(counts_sel)
cpm[1:5,1:5]
is.exprs <- rowSums(cpm>0.5) >= 3 # Apply filtering criteria
head(is.exprs)
final_counts <- counts_sel[is.exprs, ]
dim(final_counts) # 39514    20
final_counts[1:5,1:5]

# Normalize the matrix
x <- DGEList(counts=final_counts)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)

vMat <- v$E
dim(vMat) # 39514    20
vMat[1:5,1:5] 

# Keep highly variable genes:
library(matrixStats)
gene_variances <- rowVars(vMat)
head(gene_variances)

gene_variances_df <- data.frame(gene = rownames(vMat),
                               variance = gene_variances)
head(gene_variances_df)
dim(gene_variances_df) # 39514     2

top_genes <- gene_variances_df %>%
  arrange(desc(variance)) %>%
  head(5000)

top_5000_genes_mat <- vMat[rownames(vMat) %in% top_genes$gene,]
dim(top_5000_genes_mat) # 5000  20
top_5000_genes_mat[1:5,1:5]


# Transpose the matrix:
counts_t <- t(top_5000_genes_mat)
counts_t[1:5,1:5]
dim(counts_t) # 20 5000
head(metadata)

# Merge the metafile & counts:
merged_df <- merge(metadata, counts_t, by = 0)
dim(merged_df)  # 20 5004
merged_df[1:5,1:5]
merged_df <- tibble::column_to_rownames(merged_df,var = 'Row.names')
merged_df[1:4,1:6]

# Proceed with the PCA analysis:
res.pca <- PCA(merged_df[,4:dim(merged_df)[2]],  graph = FALSE)
str(merged_df)
head(colnames(merged_df))
merged_df[,c('Donor','tissue')] <- lapply(merged_df[,c('Donor','tissue')], factor)

# PCA plot for Donor:
pdf('PCA_plot_donor.pdf', width = 9, height = 7)
colors <- c(             # for the donor
  "red1",
  "#33FF57",
  "#2066a8",
  "#ff73b6",
  "#FF8C00",
  "#00BFFF",
  "#A52A2A",
  "#8ec1da",
  "yellow",
  "#c701ff"
)

p <- fviz_pca_ind(res.pca,
                  geom = "point",
                  pointshape = 19,
                  pointsize = 10,
                  alpha.ind = 0.6,
                  label = "none",
                  title = 'PCA plot for Donor',
                  habillage = merged_df$Donor,
                  palette = colors,
                  axes = c(1, 2),  # Specify the axes (PC1 and PC2)
                  addEllipses = FALSE,
                  mean.point=F)

p + theme(axis.title.x = element_text(size = 30,face = "bold"),
          axis.title.y = element_text(size = 30,face = "bold"),
          axis.text.x = element_text(size = 28,colour = "black"),
          axis.text.y = element_text(size = 28,colour = "black"),
          plot.title = element_text(size = 32,hjust = 0.5,face = "bold"),
          legend.text = element_text(size = 26,colour = "black"),
          legend.title = element_text(size = 28,face = "bold"),
          legend.key.size = unit(2.5, "lines")) + labs(color = "Donors")
dev.off()

# PCA Plot for Tissue:
pdf('PCA_plot_tissue.pdf', width = 10, height = 7)
colors2 <- c(                  # for the tissue
  "#003a7d",
  "#4ecb8d"
)

p <- fviz_pca_ind(res.pca,
                  geom = "point",
                  pointshape = 19,
                  pointsize = 10,
                  alpha.ind = 0.6,
                  label = "none",
                  title = 'PCA plot for Tissue',
                  habillage = merged_df$tissue,
                  palette = colors2,
                  axes = c(1, 2),  # Specify the axes (PC1 and PC2)
                  addEllipses = FALSE,
                  mean.point=F)

p + theme(axis.title.x = element_text(size = 30,face = "bold"),
          axis.title.y = element_text(size = 30,face = "bold"),
          axis.text.x = element_text(size = 28,colour = "black"),
          axis.text.y = element_text(size = 28,colour = "black"),
          plot.title = element_text(size = 32,hjust = 0.5,face = "bold"),
          legend.text = element_text(size = 26,colour = "black"),
          legend.title = element_text(size = 28,face = "bold"),
          legend.key.size = unit(2.5, "lines"))  + labs(color = "Tissue")
dev.off()


#----------------------------------------------------------------------------------------------#
# DIFFERENTIAL EXPRESSION ANALYSIS
#----------------------------------------------------------------------------------------------#

# Load required libraries
library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)
library(tidyverse)
library(data.table)
library(gprofiler2)
library(EnhancedVolcano)
library(edgeR)

# Load the metadata: 
metadata <- as.data.frame(read.csv("metadata.csv",header = T,row.names = 1))
metadata[1:5,]
dim(metadata) #  20  3
colnames(metadata)

table(metadata$tissue) 
# breast_tumor      normal_breast_tissue 
#       10                   10

# Read the raw counts matrix 
count <- as.data.frame(read.csv("raw_counts.csv",header = T,row.names = 1))
count[1:5,1:5]
dim(count)  # 62700    20

# Check the common samples with metadata & keep same samples as in the metadata:
length(intersect(colnames(count),rownames(metadata))) #  20  (all samples match)
setdiff(rownames(metadata),colnames(count))  # 0

count <- count[,colnames(count) %in% rownames(metadata)]
dim(count) # 62700     20
count[1:5,1:5]
all(colnames(count)==rownames(metadata)) # TRUE

# Add gene names corresponding to the ensembl ids using gprofiler2 library
all_gene_ids <- rownames(count)
head(all_gene_ids)

# Remove digits after dots
all_gene_ids_2 <- gsub("\\..*","",all_gene_ids)
head(all_gene_ids_2)

genes <- gprofiler2::gconvert(all_gene_ids_2,organism = 'hsapiens')
dim(genes) # 61481     7
colnames(genes)
gene_names=genes[,c("input","name")]  # fetching just the id and gene name columns
head(gene_names)

# Merge to add gene names in the count dataframe:
count[1:5,1:5]
count <- tibble::rownames_to_column(count,'gene_id')
count$gene_id <- gsub("\\..*","",count$gene_id)
count[1:5,1:5]

merged_counts <- merge(gene_names,count,by.x = 'input',by.y = "gene_id")
dim(merged_counts) # 61481    22
merged_counts[1:5,1:5]
table(duplicated(merged_counts$name))  # check for duplicate genes
# FALSE  TRUE 
# 41095 20386

# Remove these duplicate genes by variance:
ens_gene=merged_counts[,1:2] # assigning gene names and ensembl ids column to a new variable
head(ens_gene)

# Calculate variance of each gene
ens_gene$Variance <- apply(merged_counts[,3:dim(merged_counts)[2]], 1, var) 
dim(ens_gene) # 61481     3

head(ens_gene[duplicated(ens_gene$name),]) # check for duplicate genes

# Sample data frame
result <- ens_gene %>%
  group_by(name) %>%
  arrange(desc(Variance)) %>%
  slice(1) %>%
  ungroup()

# Print the result
dim(result) # 41095   3
head(result)

# Only keep these genes with high variance in the matrix:
final_counts=merged_counts[merged_counts$input %in% result$input,]
dim(final_counts) # 41095    22
final_counts[1:5,1:5]
table(duplicated(final_counts$name))
# FALSE 
# 41095 

# Remove ensembl id column & make gene names as rownames:
rownames(final_counts)<-NULL
final_counts <- final_counts %>% 
  tibble::column_to_rownames(var = 'name') %>%
  select(-input)
dim(final_counts) # 41095    20
final_counts[1:5,1:5]

# Now calculate cpm 
final_counts <- as.data.frame(final_counts)
cpm <- cpm(final_counts)
cpm[1:5,1:5]

# Apply the filter based on cpm:
is.exprs <- rowSums(cpm>0.5) >= 3 # Apply filtering criteria
head(is.exprs)
final_counts2 <- final_counts[is.exprs, ]
dim(final_counts2) # 28382    20
final_counts2[1:5,1:5]

# Make the order of samples same in both the files 
head(metadata)
all(rownames(metadata)==colnames(final_counts2)) # TRUE
com_exp_mat=final_counts2   # assigning to a new variable
com_metafile=metadata

dim(com_exp_mat)  # 28382   20
dim(com_metafile) #  20      3
head(com_metafile)
table(com_metafile$tissue)

group <- factor(com_metafile$tissue)
table(group)
class(group)

d <- DGEList(counts=com_exp_mat, group=group)  # group here is tissue column

# Normalize the data
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

# Create modelDesign matrix
modelDesign <- model.matrix(~0 + group)
head(modelDesign)
dim(modelDesign) # 20  2

table(modelDesign[,1])
table(modelDesign[,2])

# Make all the comparisons at this step here:
head(modelDesign)
contrast_one <- makeContrasts(
  TumorVsNormal=groupbreast_tumor-groupnormal_breast_tissue,
  levels=colnames(modelDesign))

contrast_one

fit_glm <- glmFit(d,modelDesign)
onevsrest <- glmLRT(fit_glm , contrast = contrast_one[,1])
tt_onevsrest <- topTags(onevsrest,n=nrow(d))
dim(tt_onevsrest) # 28382   5
  
# Extract the toptable:
toptable=tt_onevsrest$table
head(toptable)
dim(toptable)
write.csv(toptable,file="toptable_TumorVsNormal.csv")

# Fetch the upregulated & downregulated genes & save them:
upregulated_genes=toptable[toptable$logFC > 2 & toptable$FDR < 0.05,]
cat("Upregulated genes", ":", dim(upregulated_genes)[1], "\n")  # 370
write.csv(upregulated_genes,file="Tumor_vs_normal_upregulated_genes.csv")

downregulated_genes=toptable[toptable$logFC < -2 & toptable$FDR < 0.05,]
cat("Downregulated genes", ":", dim(downregulated_genes)[1], "\n")  # 4136
write.csv(downregulated_genes,file="Tumor_vs_normal_downregulated_genes.csv")
  
# Volcano plot:
g <- EnhancedVolcano(toptable,
                     lab = rownames(toptable),
                     selectLab = c("MIR3648-1","H3C15","GFRA1","SERPINB12","MMP9"),
                     title = "Tumor vs Normal",
                     x = 'logFC',
                     y = 'FDR',
                     drawConnectors = TRUE,
                     FCcutoff=2,
                     pCutoff = 0.05,
                     labSize = 6,
                     boxedLabels = TRUE,
                     pointSize = 3
)

ggsave(filename = "TumorVsNormal_volcano.png",
       plot = g,
       width = 8,
       height = 8,
       dpi = 500)


# Heatmap
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(circlize)

dim(com_exp_mat) # 28382    20
dim(com_metafile) # 20   3

x <- DGEList(counts=com_exp_mat)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)
vMat <- v$E
dim(vMat) # 28382    20
vMat[1:5,1:5] 
# In variance_stabilised_counts the rows are the variables and the columns correspond to the samples

# Keeping top 50 up and downregulated genes in the matrix
up <-rownames(upregulated_genes)[1:50]
down <- rownames(downregulated_genes)[1:50]
all <- c(up,down)

vMat_sel=vMat[all,]
dim(vMat_sel) # 100  20
vMat_sel[1:4,1:5]
vMat_sel=na.omit(vMat_sel)
dim(vMat_sel)  # 100  20

t_sel=t(vMat_sel)  # transpose so rows are samples and columns correspond to genes
dim(t_sel)
ti=merge(com_metafile,t_sel,by=0)
ti=ti[order(ti$tissue),,drop=F]
rownames(ti)=ti$Row.names
meta1=ti[,1:4]
exp=ti[,5:dim(ti)[2]]
dim(meta1)  # 20  4
dim(exp)  # 20  100
exp[1:4,1:4]
meta1[1:4,1:4]

pdf("DEGS_heatmap.pdf",width=8,height = 7)
Heatmap(scale(exp),name = "Scaled expression",show_row_names = TRUE, show_column_names = FALSE)
dev.off()

# For adding annotation bar
meta1=meta1[order(meta1$tissue),]
exp=exp[rownames(meta1),]
head(meta1)

# Assign colours
table(meta1$tissue)
color_tissues <- c(                  # for the tissue
  "breast_tumor"="#003a7d",
  "normal_breast_tissue"= "#4ecb8d"
)

meta1$tissue <- factor(meta1$tissue, levels = c("breast_tumor","normal_breast_tissue"))

pdf("DEGS_heatmap_annotated.pdf",width=8,height = 7)
h1=Heatmap(meta1$tissue,cluster_rows = F,width = unit(1, "cm"),name="Tissue",show_row_names = F,col=color_tissues)
h2=Heatmap(scale(exp),name="Scaled Expression",show_column_names = FALSE,
           row_names_gp = grid::gpar(fontsize = 8),cluster_columns = FALSE)
h1+h2
dev.off()


#----------------------------------------------------------------------------------------------#
# OVER-REPRESENTATION ANALYSIS (ORA)
#----------------------------------------------------------------------------------------------#

# Read the hallmark gmt file:
gmt_file <- upload_GMT_file("hallmark.gmt") 
gmt_file[[1]] 

# Define the background genes:
background_genes <- rownames(com_exp_mat)
length(background_genes) # 28382
head(background_genes)
write.csv(background_genes,"background_genes_28382.csv")

# For Upregulated Genes:
gostres_up <- gost(query = rownames(upregulated_genes), organism = gmt_file[[1]],
                     ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
                     exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE,
                     user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated",
                     custom_bg = background_genes, numeric_ns = "", sources = NULL, as_short_link = FALSE)
  
result_up <- as.data.frame(gostres_up$result)
result_up$query <- "Upregulated"

# Printing the most significant pathway:
most_significant_pathway <- result_up[which.min(result_up$p_value), ]$term_id
cat("Most significant pathway", ":", most_significant_pathway, "\n")

  
# For Downregulated Genes:
gostres_down <- gost(query = rownames(downregulated_genes), organism = gmt_file[[1]],
                     ordered_query = FALSE, multi_query = FALSE, significant = TRUE,
                     exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE,
                     user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated",
                     custom_bg = background_genes, numeric_ns = "", sources = NULL, as_short_link = FALSE)

result_down <- as.data.frame(gostres_down$result)
result_down$query <- "Downregulated"
  
# Printing the most significant pathway:
most_significant_pathway <- result_down[which.min(result_down$p_value), ]$term_id
cat("Most significant pathway",":", most_significant_pathway, "\n")

# Save the final results:
result_up2 <- apply(result_up,2,as.character)
write.csv(result_up2, file = "Pathway_results_upregulated.csv", row.names = FALSE)

result_down2 <- apply(result_down,2,as.character)
write.csv(result_down2, file = "Pathway_results_downregulated.csv", row.names = FALSE)

# Plot
# Combine results
all_results <- rbind(result_up, result_down)
dim(all_results) # 13  16
all_results$term_id <- gsub("HALLMARK_","",all_results$term_id)

# Keep top 3 pathways per query based on p-value
top_results <- all_results %>%
  group_by(query) %>%
  slice_min(order_by = p_value, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    query = factor(query, levels = c("Upregulated","Downregulated")),
    log10_p = -log10(p_value)
  ) %>%
  arrange(query, desc(log10_p)) %>%
  mutate(term_id = factor(term_id, levels = rev(unique(term_id))))

ggplot(top_results, aes(x = term_id, y = log10_p, fill = query)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Pathway",
    y = "-log10(p-value)",
    title = "Top 3 Enriched Pathways per Query"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Upregulated" = "#E41A1C",
                               "Downregulated" = "#377EB8")) +
  theme(axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 16))


#----------------------------------------------------------------------------------------------#
# Gene Set Variation Analysis (GSVA)
#----------------------------------------------------------------------------------------------#

library(qusage)
library(GSVA)
library(ActivePathways)
library(ComplexHeatmap)
library(circlize)

# Continue with the files we prepared above
dim(com_exp_mat)  # 28382   20
dim(com_metafile) #  20     3

# Normalize the matrix
x <- DGEList(counts=com_exp_mat)
x <- calcNormFactors(x, method='TMM')
v <- voom(x, plot=T)

# Extract normalized matrix
vMat <- v$E
dim(vMat) # 28382    20
vMat[1:5,1:5]  # feed this matrix into gsva 

# Read the hallmark gmt file
gmt<- read.gmt("hallmark.gmt")

# Run GSVA
d <- gsva(vMat,gmt, method = "gsva",kcdf= "Gaussian", min.sz= 2, verbose = T)
class(d)
d <- as.data.frame(d)

# Save the results
write.csv(d,file="GSVA_result.csv")

# For visualizing the results through a heatmap
gsva_scores=t(d)
head(gsva_scores)
com_metafile[1:4,]
finaldf=merge(com_metafile,gsva_scores,by=0)
all(com_metafile$Sample==finaldf$Sample) # TRUE

str(finaldf)
finaldf$tissue=as.factor(finaldf$tissue)
finaldf[1:5,1:4]
finaldf <- tibble::column_to_rownames(finaldf,var = 'Row.names')
dim(finaldf) #  20  53

ti=finaldf[order(finaldf$tissue),,drop=F]
meta1=ti[,1:3]
exp=ti[,4:dim(ti)[2]]
dim(meta1)  # 20  3
dim(exp)    # 20  50
exp[1:4,1:4]
head(meta1)
exp <- as.matrix(exp)
dim(exp)  # 20  50
colnames(exp)<-gsub("HALLMARK_","",colnames(exp))
exp[1:4,1:4]

pdf("gsva_scores_heatmap.pdf",width=16,height = 15)
ht <- Heatmap(t(scale(exp)),cluster_rows = T,cluster_columns = T,name="Scaled Expression",show_column_names =F)
draw(ht,padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

## For adding annotation:
# meta1=meta1[order(meta1$tissue),]
# exp=exp[rownames(meta1),]
# head(meta1)
# 
# # Assign colours
# table(meta1$tissue)
# color_tissues <- c(                  # for the tissue
#   "breast_tumor"="#003a7d",
#   "normal_breast_tissue"= "#4ecb8d"
# )
# 
# meta1$tissue <- factor(meta1$tissue, levels = c("breast_tumor","normal_breast_tissue"))
# 
# h1=Heatmap(meta1$tissue,cluster_rows = F,width = unit(1, "cm"),name="Tissue",
#            column_names_gp = grid::gpar(fontsize=22,fontface="bold"),
#            show_row_names = F,
#            col=color_tissues,
#            heatmap_legend_param = list(
#              labels_gp = gpar(fontsize = 18),
#              title_gp = gpar(fontsize = 18,fontface = 'bold')))
# 
# h2=Heatmap(scale(exp),name="Scaled Expression",column_names_gp = grid::gpar(fontsize = 15,fontface = 'bold'),
#            row_names_gp = grid::gpar(fontsize = 15,fontface = 'bold'),
#            heatmap_legend_param = list(labels_gp = gpar(fontsize = 18),
#                                                       title_gp = gpar(fontsize=18,fontface = 'bold')))
# 
# pdf("gsva_scores_heatmap.pdf",width=12,height = 15)
# # h1+h2
# draw(
#   h1 + h2,
#   padding = unit(c(10, 10, 10, 10), "mm")
# )
# dev.off()


#----------------------------------------------------------------------------------------------#
# Gene Set Enrichment Analysis (GSEA)
#----------------------------------------------------------------------------------------------#

# BiocManager::install("fgsea")
library(fgsea)
library(dplyr)
library(tidyverse)
library(qusage)
library(writexl)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(data.table)

# Load the hallmark gmt file & convert to a binary matrix:
gmt_file <- gmtPathways('hallmark.gmt') 
names(gmt_file) 
total_unique_genes <- unique(unlist(gmt_file))
length(total_unique_genes)  # 4384

# Convert gmt file to a matrix with the genes as rows and for each pathway (columns) the values are 0 or 1
mat <- matrix(NA, dimnames = list(total_unique_genes, names(gmt_file)),
              nrow = length(total_unique_genes), ncol = length(gmt_file))

for (i in 1:dim(mat)[2]){
  mat[,i] <- as.numeric(total_unique_genes %in% gmt_file[[i]])
}

# Read the genes given for DEG analysis which will be used as background genes here:
background_genes <- as.data.frame(read.csv('background_genes_28382.csv',row.names = 1))
dim(background_genes) #  28382    1
head(background_genes)
background_genes <- background_genes$x

# Subset to the genes that are present in our data to avoid bias:
com <- intersect(background_genes, total_unique_genes)
length(com) # 4278
mat <- mat[com, colnames(mat)[which(colSums(mat[com,])>5)]] # filter for gene sets with more than 5 genes annotated
dim(mat) # 4278   50

# Define a function to convert this matrix back to a list:
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}


# And get the list again using the function we previously defined
final_pathway_list <- matrix_to_list(mat)

# Now prepare ranked list of DEGs
upreg <- read.csv('Tumor_vs_normal_upregulated_genes.csv',header = T,row.names = 1)
downreg <- read.csv('Tumor_vs_normal_downregulated_genes.csv',header = T,row.names = 1)
deg_all <- rbind(upreg,downreg)
head(deg_all)
deg_all<-as.data.frame(deg_all)
deg_all<-tibble::rownames_to_column(deg_all,"gene")
head(deg_all)

# deg_all$FDR[deg_all$FDR==0] <- 1e-300
rankings <- sign(deg_all$logFC)*(-log10(deg_all$FDR)) # we will use both logfc and pvalue
names(rankings) <- deg_all$gene # gene names chosen here
head(rankings)
#    MMP11     H3C15     PRAME MIR3648-2 MIR3648-1 SMYD3-IT1 
# 6.425976  5.253664  4.864155  4.669878  4.499926  4.388352 

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)
max(rankings)
min(rankings)

# Run GSEA
GSEAres <- fgsea(pathways = final_pathway_list, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std',
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

# Check results
head(GSEAres[order(pval), ])
sum(GSEAres[, padj < 0.05])

# Save the results
str(GSEAres)
GSEAres2 <- as.data.frame(apply(GSEAres,2,as.character))
write_xlsx(GSEAres2,'Tumor_Vs_Normal_GSEA.xlsx')

# Formatting:
GSEAres_converted <- GSEAres %>%
  mutate(
    pval = as.numeric(pval),
    padj = as.numeric(padj),
    log2err = as.numeric(log2err),
    ES = as.numeric(ES),
    NES = as.numeric(NES),
    size = as.integer(size),
    leadingEdge = sapply(leadingEdge, paste, collapse = ", ")  # Convert leadingEdge back to a list of character vectors
  )

# Convert to data.table
setDT(GSEAres_converted)

# Check the structure of the new data frame or data.table
str(GSEAres_converted)

# Extract top pathways:
topPathwaysUp <- GSEAres_converted[NES > 0][head(order(padj)), pathway]
topPathwaysDown <- GSEAres_converted[NES < 0][head(order(padj)),pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
topPathways

# Select only independent pathways, removing redundancies/similar pathways
collapsedPathways <- collapsePathways(GSEAres_converted[order(GSEAres_converted)][padj < 0.05], final_pathway_list, rankings)
mainPathways <- GSEAres_converted[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

pdf('tumor_vs_normal_GSEA_table.pdf',width = 12,height = 8)
plotGseaTable(final_pathway_list[topPathways],
              stats = rankings, 
              fgseaRes = GSEAres_converted, 
              gseaParam = 0.5)
dev.off()

# Plot the most significantly enriched pathway
plotEnrichment(final_pathway_list[[head(GSEAres_converted[order(padj), ], 1)$pathway]],
               rankings,ticksSize = 0.5) +
  labs(title = head(GSEAres_converted[order(padj), ], 1)$pathway)

# Plot a specific pathway
plotEnrichment(final_pathway_list[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
               rankings, ticksSize = 0.5) + 
  labs(title = "EPITHELIAL MESENCHYMAL TRANSITION") +
  xlab(label = "Rank") +
  ylab(label = "Enrichment Score") +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5,face = "bold",family = "Arial",colour = "black"),  
    axis.title.x = element_text(size = 9, face = "bold",family = "Arial",colour = "black"),  
    axis.title.y = element_text(size = 9, face = "bold",family = "Arial",colour = "black"),  
    axis.text.x = element_text(size = 9, face = "bold",family = "Arial",colour = "black"),
    axis.text.y = element_text(size = 9, face = "bold",family = "Arial",colour = "black"))

# For making a barplot representation of results
# Just fetch the top pathways:
length(topPathways) # 12
GSEAres2 <- GSEAres_converted[GSEAres_converted$pathway %in% topPathways,]
dim(GSEAres2) # 12  8

GSEAres2$pathway <- gsub("HALLMARK_","",GSEAres2$pathway)
GSEAres2$neglog10adjpval <- -log10(GSEAres2$padj)
GSEAres2$color <- ifelse(GSEAres2$NES > 0, "#ff9d3a", "#00b0be")

# Create a bar plot
g<-ggplot(GSEAres2, aes(x = reorder(pathway, NES), y = NES, fill = color)) +
  geom_bar(stat = "identity") +  
  scale_fill_identity() +         
  coord_flip() +                  
  theme_gray() +              
  labs(title = "GSEA Enrichment Plot", 
       x = "Pathway", 
       y = "Normalized Enrichment Score (NES)") +
  theme(axis.text.x = element_text(size = 18, hjust = 1,colour = "black",face = "bold"),
        plot.title = element_text(size = 24,face = "bold",colour = "black",hjust = 0.5),
        axis.title.x = element_text(size = 20,face = "bold",colour = "black"),
        axis.title.y = element_text(size = 20,face = "bold",colour = "black"),
        axis.text.y = element_text(size = 18,colour = "black",face = "bold")) +
  geom_text(aes(label = sprintf("%.2f", neglog10adjpval)),  
            position = position_stack(vjust = 0.5), 
            size = 9, color = "black",fontface="bold") 
g
ggsave(filename = "gsea_barplot_tumor_vs_normal.png", g, width = 12, height = 9, dpi = 500)

############################################################################################
############################################################################################
