## Library

library(cellrangerRkit)
library(devtools)
library(Seurat)
library(Matrix)
library(monocle)
library(DDRTree)
library(pheatmap)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)
library(clues)

## Setup

# load data
pipestance_path1 <- "../ED17/outs/"
med165 <- load_cellranger_matrix(pipestance_path1)
pipestance_path2 <- "../PND7/outs/"
mp7 <- load_cellranger_matrix(pipestance_path2)
pipestance_path3 <- "../PND14/outs/"
mp14 <- load_cellranger_matrix(pipestance_path3)

# rename gene symbol column
my_feat16 <- fData(med165)
names(my_feat16) <- c('id', 'gene_short_name')
my_feat7 <- fData(mp7)
names(my_feat16) <- c('id', 'gene_short_name')
my_feat14 <- fData(mp14)
names(my_feat16) <- c('id', 'gene_short_name')

# New data set
my_cds16 <- newCellDataSet(exprs(med165),
                           phenoData = new("AnnotatedDataFrame", data = pData(med165)),
                           featureData = new("AnnotatedDataFrame", data = my_feat16),
                           lowerDetectionLimit = 0.5,
                           expressionFamily = negbinomial.size())
my_cds16 <- estimateSizeFactors(my_cds16)
my_cds16 <- estimateDispersions(my_cds16)
my_cds16 <- detectGenes(my_cds16, min_expr = 0.1)

my_cds7 <- newCellDataSet(exprs(mp7),
                           phenoData = new("AnnotatedDataFrame", data = pData(mp7)),
                           featureData = new("AnnotatedDataFrame", data = my_feat7),
                           lowerDetectionLimit = 0.5,
                           expressionFamily = negbinomial.size())
my_cds7 <- estimateSizeFactors(my_cds7)
my_cds7 <- estimateDispersions(my_cds7)
my_cds7 <- detectGenes(my_cds7, min_expr = 0.1)

my_cds14 <- newCellDataSet(exprs(mp14),
                           phenoData = new("AnnotatedDataFrame", data = pData(mp14)),
                           featureData = new("AnnotatedDataFrame", data = my_feat14),
                           lowerDetectionLimit = 0.5,
                           expressionFamily = negbinomial.size())
my_cds14 <- estimateSizeFactors(my_cds14)
my_cds14 <- estimateDispersions(my_cds14)
my_cds14 <- detectGenes(my_cds14, min_expr = 0.1)

# standardize to Z-distribution
x16 <- pData(my_cds16)$num_genes_expressed
x16_1 <- (x16 - mean(x16)) / sd(x16)
summary(x16_1)
df16 <- data.frame(x16 = x16_1)
ggplot(df16, aes(x16)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')

xp7 <- pData(my_cds7)$num_genes_expressed
xp7_1 <- (xp7 - mean(xp7)) / sd(xp7)
summary(xp7_1)
dfp7 <- data.frame(xp7 = xp7_1)
ggplot(dfp7, aes(xp7)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')

xp14 <- pData(my_cds14)$num_genes_expressed
xp14_1 <- (xp14 - mean(xp14)) / sd(xp14)
summary(xp14_1)
dfp14 <- data.frame(xp14 = xp14_1)
ggplot(dfp14, aes(xp14)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')

# add a UMI column into phenoData
pData(my_cds16)$UMI <- Matrix::colSums(exprs(my_cds16))
ggplot(pData(my_cds16), aes(num_genes_expressed, UMI)) + geom_point()
disp_table16 <- dispersionTable(my_cds16)
table(disp_table16$mean_expression>=0.1)
unsup_clustering_genes16 <- subset(disp_table16, mean_expression >= 0.1)
my_cds16 <- setOrderingFilter(my_cds16, unsup_clustering_genes16$gene_id)
plot_ordering_genes(my_cds16)
plot_pc_variance_explained(my_cds16, return_all = FALSE)

pData(my_cds7)$UMI <- Matrix::colSums(exprs(my_cds7))
ggplot(pData(my_cds7), aes(num_genes_expressed, UMI)) + geom_point()
disp_tablep7 <- dispersionTable(my_cds7)
table(disp_tablep7$mean_expression>=0.1)
unsup_clustering_genesp7 <- subset(disp_tablep7, mean_expression >= 0.1)
my_cds7 <- setOrderingFilter(my_cds7, unsup_clustering_genesp7$gene_id)
plot_ordering_genes(my_cds7)
plot_pc_variance_explained(my_cds7, return_all = FALSE)

pData(my_cds14)$UMI <- Matrix::colSums(exprs(my_cds14))
ggplot(pData(my_cds14), aes(num_genes_expressed, UMI)) + geom_point()
disp_tablep14 <- dispersionTable(my_cds14)
table(disp_tablep14$mean_expression>=0.1)
unsup_clustering_genesp14 <- subset(disp_tablep14, mean_expression >= 0.1)
my_cds14 <- setOrderingFilter(my_cds14, unsup_clustering_genesp14$gene_id)
plot_ordering_genes(my_cds14)
plot_pc_variance_explained(my_cds14, return_all = FALSE)

## Clustering
# reduceDimension, num_dim = 10
my_cds16 <- reduceDimension(my_cds16, max_components = 2, num_dim = 10,
                            reduction_method = 'tSNE', verbose = TRUE)
my_cds16 <- clusterCells(my_cds16, num_clusters = 15)
my_cluster_dim_10 <- pData(my_cds16)$Cluster
plot_cell_clusters(my_cds16)

my_cds7 <- reduceDimension(my_cds7, max_components = 2, num_dim = 10,
                            reduction_method = 'tSNE', verbose = TRUE)
my_cds7 <- clusterCells(my_cds7, num_clusters = 15)
my_cluster_dim_10_p7 <- pData(my_cds7)$Cluster
plot_cell_clusters(my_cds7)

my_cds14 <- reduceDimension(my_cds14, max_components = 2, num_dim = 10,
                             reduction_method = 'tSNE', verbose = TRUE)
my_cds14 <- clusterCells(my_cds14, num_clusters = 15)
my_cluster_dim_10_p14 <- pData(my_cds14)$Cluster
plot_cell_clusters(my_cds14)

# vector to phenoData
my_vector16 <- rep('no', nrow(pData(my_cds16)))
my_vector16[pData(my_cds16)$Cluster == 1] <- rep('yes', sum(pData(my_cds16)$Cluster == 1))
pData(my_cds16)$test <- my_vector16
head(pData(my_cds16))
length(unsup_clustering_genes16$gene_id)
de_cluster_one16 <- differentialGeneTest(my_cds16[unsup_clustering_genes16$gene_id,],
                                         fullModelFormulaStr = '~test',
                                         cores = 8)

my_vectorp7 <- rep('no', nrow(pData(my_cds7)))
my_vectorp7[pData(my_cds7)$Cluster == 1] <- rep('yes', sum(pData(my_cds7)$Cluster == 1))
pData(my_cds7)$test <- my_vectorp7
head(pData(my_cds7))
length(unsup_clustering_genesp7$gene_id)
de_cluster_onep7 <- differentialGeneTest(my_cds7[unsup_clustering_genesp7$gene_id,],
                                         fullModelFormulaStr = '~test',
                                         cores = 8)

my_vectorp14 <- rep('no', nrow(pData(my_cds14)))
my_vectorp14[pData(my_cds14)$Cluster == 1] <- rep('yes', sum(pData(my_cds14)$Cluster == 1))
pData(my_cds14)$test <- my_vectorp14
head(pData(my_cds14))
length(unsup_clustering_genesp14$gene_id)
de_cluster_onep14 <- differentialGeneTest(my_cds14[unsup_clustering_genesp14$gene_id,],
                                          fullModelFormulaStr = '~test',
                                          cores = 8)

# order by q-value
de_cluster_one16 %>% arrange(qval) %>% head()
de_cluster_onep7 %>% arrange(qval) %>% head()
de_cluster_onep14 %>% arrange(qval) %>% head()

# Plots - Non-somatic cell detection
plot_cell_clusters(my_cds16, markers="Dazl")
plot_cell_clusters(my_cds16, markers="Sycp1")
plot_cell_clusters(my_cds16, markers="Sycp3")
plot_cell_clusters(my_cds16, markers="Ddx4")
plot_cell_clusters(my_cds16, markers="Alas2")
plot_cell_clusters(my_cds16, markers="Hbb-bs")
plot_cell_clusters(my_cds16, markers="Lyz2")
plot_cell_clusters(my_cds16, markers="Ccl24")
plot_cell_clusters(my_cds16, markers="Cd34")
plot_cell_clusters(my_cds16, markers="Cd36")
plot_cell_clusters(my_cds16, markers="Icam2")
plot_cell_clusters(my_cds16, markers="Pecam1")
plot_cell_clusters(my_cds7, markers="Dazl")
plot_cell_clusters(my_cds7, markers="Sycp1")
plot_cell_clusters(my_cds7, markers="Sycp3")
plot_cell_clusters(my_cds7, markers="Ddx4")
plot_cell_clusters(my_cds7, markers="Alas2")
plot_cell_clusters(my_cds7, markers="Hbb-bs")
plot_cell_clusters(my_cds7, markers="Lyz2")
plot_cell_clusters(my_cds7, markers="Ccl24")
plot_cell_clusters(my_cds7, markers="Cd34")
plot_cell_clusters(my_cds7, markers="Cd36")
plot_cell_clusters(my_cds7, markers="Icam2")
plot_cell_clusters(my_cds7, markers="Pecam1")
plot_cell_clusters(my_cds14, markers="Dazl")
plot_cell_clusters(my_cds14, markers="Sycp1")
plot_cell_clusters(my_cds14, markers="Sycp3")
plot_cell_clusters(my_cds14, markers="Ddx4")
plot_cell_clusters(my_cds14, markers="Alas2")
plot_cell_clusters(my_cds14, markers="Hbb-bs")
plot_cell_clusters(my_cds14, markers="Lyz2")
plot_cell_clusters(my_cds14, markers="Ccl24")
plot_cell_clusters(my_cds14, markers="Cd34")
plot_cell_clusters(my_cds14, markers="Cd36")
plot_cell_clusters(my_cds14, markers="Icam2")
plot_cell_clusters(my_cds14, markers="Pecam1")

# Plots - Cell types
plot_cell_clusters(my_cds16)
plot_cell_clusters(my_cds16, markers="Lgr5")
plot_cell_clusters(my_cds16, markers="Esr1")
plot_cell_clusters(my_cds16, markers="Esr2")
plot_cell_clusters(my_cds16, markers="Wt1")
plot_cell_clusters(my_cds16, markers="Foxl2")
plot_cell_clusters(my_cds16, markers="Nr2f2")
plot_cell_clusters(my_cds16, markers="Bmp15")
plot_cell_clusters(my_cds16, markers="Gdf9")
plot_cell_clusters(my_cds16, markers="Tgfb1")
plot_cell_clusters(my_cds16, markers="Tgfbr1")
plot_cell_clusters(my_cds16, markers="Tgfbr2")
plot_cell_clusters(my_cds16, markers="Tgfbr3")
plot_cell_clusters(my_cds16, markers="Inhba")
plot_cell_clusters(my_cds16, markers="Kit")
plot_cell_clusters(my_cds16, markers="Dazl")
plot_cell_clusters(my_cds16, markers="Sycp1")
plot_cell_clusters(my_cds16, markers="Sycp3")
plot_cell_clusters(my_cds16, markers="Ddx4")
plot_cell_clusters(my_cds16, markers="Gli1")
plot_cell_clusters(my_cds16, markers="Cyp11a1")
plot_cell_clusters(my_cds16, markers="Cyp17a1")
plot_cell_clusters(my_cds16, markers="Gli1")
# OSE
plot_cell_clusters(my_cds16, markers="Lgals7")
plot_cell_clusters(my_cds16, markers="Krt18")
plot_cell_clusters(my_cds16, markers="Krt8")
plot_cell_clusters(my_cds16, markers="Ly6e")
# GC
plot_cell_clusters(my_cds16, markers="Fst")
plot_cell_clusters(my_cds16, markers="Amh")
plot_cell_clusters(my_cds16, markers="Inha")
plot_cell_clusters(my_cds16, markers="Ube2c")
plot_cell_clusters(my_cds16, markers="Car14")
plot_cell_clusters(my_cds16, markers="Kctd14")
# Theca
plot_cell_clusters(my_cds16, markers="Hsd3b1")
plot_cell_clusters(my_cds16, markers="Col1a2")
plot_cell_clusters(my_cds16, markers="Loxl2")
plot_cell_clusters(my_cds16, markers="Ramp2")
plot_cell_clusters(my_cds16, markers="Ogn")
plot_cell_clusters(my_cds16, markers="Dcn")

plot_cell_clusters(my_cds7)
plot_cell_clusters(my_cds7, markers="Lgr5")
plot_cell_clusters(my_cds7, markers="Esr1")
plot_cell_clusters(my_cds7, markers="Esr2")
plot_cell_clusters(my_cds7, markers="Wt1")
plot_cell_clusters(my_cds7, markers="Foxl2")
plot_cell_clusters(my_cds7, markers="Nr2f2")
plot_cell_clusters(my_cds7, markers="Bmp15")
plot_cell_clusters(my_cds7, markers="Gdf9")
plot_cell_clusters(my_cds7, markers="Tgfb1")
plot_cell_clusters(my_cds7, markers="Tgfbr1")
plot_cell_clusters(my_cds7, markers="Tgfbr2")
plot_cell_clusters(my_cds7, markers="Tgfbr3")
plot_cell_clusters(my_cds7, markers="Inhba")
plot_cell_clusters(my_cds7, markers="Kit")
plot_cell_clusters(my_cds7, markers="Dazl")
plot_cell_clusters(my_cds7, markers="Sycp1")
plot_cell_clusters(my_cds7, markers="Sycp3")
plot_cell_clusters(my_cds7, markers="Ddx4")
plot_cell_clusters(my_cds7, markers="Gli1")
plot_cell_clusters(my_cds7, markers="Cyp11a1")
plot_cell_clusters(my_cds7, markers="Cyp17a1")
plot_cell_clusters(my_cds7, markers="Gli1")
# OSE
plot_cell_clusters(my_cds7, markers="Lgals7")
plot_cell_clusters(my_cds7, markers="Krt18")
plot_cell_clusters(my_cds7, markers="Krt8")
plot_cell_clusters(my_cds7, markers="Ly6e")
# GC
plot_cell_clusters(my_cds7, markers="Fst")
plot_cell_clusters(my_cds7, markers="Amh")
plot_cell_clusters(my_cds7, markers="Inha")
plot_cell_clusters(my_cds7, markers="Ube2c")
plot_cell_clusters(my_cds7, markers="Car14")
plot_cell_clusters(my_cds7, markers="Kctd14")
# Theca
plot_cell_clusters(my_cds7, markers="Hsd3b1")
plot_cell_clusters(my_cds7, markers="Col1a2")
plot_cell_clusters(my_cds7, markers="Loxl2")
plot_cell_clusters(my_cds7, markers="Ramp2")
plot_cell_clusters(my_cds7, markers="Ogn")
plot_cell_clusters(my_cds7, markers="Dcn")

plot_cell_clusters(my_cds14)
plot_cell_clusters(my_cds14, markers="Lgr5")
plot_cell_clusters(my_cds14, markers="Esr1")
plot_cell_clusters(my_cds14, markers="Esr2")
plot_cell_clusters(my_cds14, markers="Wt1")
plot_cell_clusters(my_cds14, markers="Foxl2")
plot_cell_clusters(my_cds14, markers="Nr2f2")
plot_cell_clusters(my_cds14, markers="Bmp15")
plot_cell_clusters(my_cds14, markers="Gdf9")
plot_cell_clusters(my_cds14, markers="Tgfb1")
plot_cell_clusters(my_cds14, markers="Tgfbr1")
plot_cell_clusters(my_cds14, markers="Tgfbr2")
plot_cell_clusters(my_cds14, markers="Tgfbr3")
plot_cell_clusters(my_cds14, markers="Inhba")
plot_cell_clusters(my_cds14, markers="Kit")
plot_cell_clusters(my_cds14, markers="Dazl")
plot_cell_clusters(my_cds14, markers="Sycp1")
plot_cell_clusters(my_cds14, markers="Sycp3")
plot_cell_clusters(my_cds14, markers="Ddx4")
plot_cell_clusters(my_cds14, markers="Gli1")
plot_cell_clusters(my_cds14, markers="Cyp11a1")
plot_cell_clusters(my_cds14, markers="Cyp17a1")
plot_cell_clusters(my_cds14, markers="Gli1")
# OSE
plot_cell_clusters(my_cds14, markers="Lgals7")
plot_cell_clusters(my_cds14, markers="Krt18")
plot_cell_clusters(my_cds14, markers="Krt8")
plot_cell_clusters(my_cds14, markers="Ly6e")
# GC
plot_cell_clusters(my_cds14, markers="Fst")
plot_cell_clusters(my_cds14, markers="Amh")
plot_cell_clusters(my_cds14, markers="Inha")
plot_cell_clusters(my_cds14, markers="Ube2c")
plot_cell_clusters(my_cds14, markers="Car14")
plot_cell_clusters(my_cds14, markers="Kctd14")
# Theca
plot_cell_clusters(my_cds14, markers="Hsd3b1")
plot_cell_clusters(my_cds14, markers="Col1a2")
plot_cell_clusters(my_cds14, markers="Loxl2")
plot_cell_clusters(my_cds14, markers="Ramp2")
plot_cell_clusters(my_cds14, markers="Ogn")
plot_cell_clusters(my_cds14, markers="Dcn")

# Isolation of ovarian somatic cells 
pData(my_cds16)$my_colour16 <- pData(my_cds16)$Cluster == 12 | pData(my_cds16)$Cluster == 14 | 
  pData(my_cds16)$Cluster == 15 | pData(my_cds16)$Cluster == 6 | 
  pData(my_cds16)$Cluster == 5 | pData(my_cds16)$Cluster == 3 | 
  pData(my_cds16)$Cluster == 9 
plot_cell_clusters(my_cds16, color_by = 'my_colour16')
pData(my_cds7)$my_colour7 <- pData(my_cds7)$Cluster == 1 | pData(my_cds7)$Cluster == 2 | 
  pData(my_cds7)$Cluster == 3 | pData(my_cds7)$Cluster == 4 | pData(my_cds7)$Cluster == 6 | 
  pData(my_cds7)$Cluster == 8 | pData(my_cds7)$Cluster == 9 | pData(my_cds7)$Cluster == 10 | 
  pData(my_cds7)$Cluster == 11 | pData(my_cds7)$Cluster == 13 | pData(my_cds7)$Cluster == 14 | 
  pData(my_cds7)$Cluster == 15 
plot_cell_clusters(my_cds7, color_by = 'my_colour7')
pData(my_cds14)$my_colour14 <- pData(my_cds14)$Cluster == 1 | pData(my_cds14)$Cluster == 2 | 
  pData(my_cds14)$Cluster == 3 | pData(my_cds14)$Cluster == 4 | pData(my_cds14)$Cluster == 5 | 
  pData(my_cds14)$Cluster == 6 | pData(my_cds14)$Cluster == 7 | pData(my_cds14)$Cluster == 8 | 
  pData(my_cds14)$Cluster == 9 | pData(my_cds14)$Cluster == 10 | pData(my_cds14)$Cluster == 12 | 
  pData(my_cds14)$Cluster == 14 
plot_cell_clusters(my_cds14, color_by = 'my_colour14')

## E16.5, PND7, PND14 combined

# load data (combined)
pipestance_pathcom <- "../Combined/outs/"
gbmcomb <- load_cellranger_matrix(pipestance_pathcom)

# rename gene symbol column
my_featcomb <- fData(gbmcomb)
names(my_featcomb) <- c('id', 'gene_short_name')

# New data set
my_cdscomb <- newCellDataSet(exprs(gbmcomb),
                             phenoData = new("AnnotatedDataFrame", data = pData(gbmcomb)),
                             featureData = new("AnnotatedDataFrame", data = my_featcomb),
                             lowerDetectionLimit = 0.5,
                             expressionFamily = negbinomial.size())
my_cdscomb <- estimateSizeFactors(my_cdscomb)
my_cdscomb <- estimateDispersions(my_cdscomb)
my_cdscomb <- detectGenes(my_cdscomb, min_expr = 0.1)

summary(fData(my_cdscomb)$num_cells_expressed)
summary(pData(my_cdscomb)$num_genes_expressed)

# standardize to Z-distribution
xcomb <- pData(my_cdscomb)$num_genes_expressed
xcomb_1 <- (xcomb - mean(xcomb)) / sd(xcomb)
summary(xcomb_1)
dfcomb <- data.frame(xcomb = xcomb_1)
ggplot(dfcomb, aes(xcomb)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')

# add a UMI column into phenoData
pData(my_cdscomb)$UMI <- Matrix::colSums(exprs(my_cdscomb))
ggplot(pData(my_cdscomb), aes(num_genes_expressed, UMI)) + geom_point()
disp_tablecomb <- dispersionTable(my_cdscomb)
table(disp_tablecomb$mean_expression>=0.1)
unsup_clustering_genescomb <- subset(disp_tablecomb, mean_expression >= 0.1)
my_cdscomb <- setOrderingFilter(my_cdscomb, unsup_clustering_genescomb$gene_id)
plot_ordering_genes(my_cdscomb)
plot_pc_variance_explained(my_cdscomb, return_all = FALSE)

# reduce Dimension, include 10 dimensions
my_cdscomb <- reduceDimension(my_cdscomb, max_components = 2, num_dim = 10,
                              reduction_method = 'tSNE', verbose = TRUE)
my_cdscomb <- clusterCells(my_cdscomb, num_clusters = 20)
plot_cell_clusters(my_cdscomb, color_by = 'Cluster')

# vector to phenoData
my_vectorcomb <- rep('no', nrow(pData(my_cdscomb)))
my_vectorcomb[pData(my_cdscomb)$Cluster == 1] <- rep('yes', sum(pData(my_cdscomb)$Cluster == 1))
pData(my_cdscomb)$test <- my_vectorcomb
length(unsup_clustering_genescomb$gene_id)
de_cluster_onecomb <- differentialGeneTest(my_cdscomb[unsup_clustering_genescomb$gene_id,],
                                           fullModelFormulaStr = '~test',
                                           cores = 8)

# order by q-value
de_cluster_onecomb %>% arrange(qval) %>% head()
plot_genes_jitter(my_cdscomb['ENSMUSG00000026238.14',], grouping = "Cluster")
plot_cell_clusters(my_cdscomb)
plot_cell_clusters(my_cdscomb,
                   color_by="Cluster",
                   cell_size = 1,
                   show_group_id = T)

# Plots - Cell types
plot_cell_clusters(my_cdscomb)
plot_cell_clusters(my_cdscomb, markers="Lgr5")
plot_cell_clusters(my_cdscomb, markers="Krt18")
plot_cell_clusters(my_cdscomb, markers="Wt1")
plot_cell_clusters(my_cdscomb, markers="Fst")
plot_cell_clusters(my_cdscomb, markers="Inha")
plot_cell_clusters(my_cdscomb, markers="Amh")
plot_cell_clusters(my_cdscomb, markers="Cyp17a1")
plot_cell_clusters(my_cdscomb, markers="Sycp1")
plot_cell_clusters(my_cdscomb, markers="Icam2")
plot_cell_clusters(my_cdscomb, markers="Alas2")
plot_cell_clusters(my_cdscomb, markers="Lyz2")
plot_cell_clusters(my_cdscomb, markers="Nr2f2")

# clusters GC lineage
pData(my_cdscomb)$my_colourcom <- pData(my_cdscomb)$Cluster == 2 | pData(my_cdscomb)$Cluster == 7 | pData(my_cdscomb)$Cluster == 13
# where TRUE means cluster membership in clusters_GC-related
# pData(my_cdscomb)$my_colourcom <- pData(my_cdscomb)$Cluster == 19 | pData(my_cdscomb)$Cluster == 4 | pData(my_cdscomb)$Cluster == 10 | pData(my_cdscomb)$Cluster == 13
# where TRUE means cluster membership in clusters_GC-related
# pData(my_cdscomb)$my_colourcom <- pData(my_cdscomb)$Cluster == 6 | pData(my_cdscomb)$Cluster == 3 | pData(my_cdscomb)$Cluster == 12 
plot_cell_clusters(my_cdscomb, color_by = 'Cluster')
plot_cell_clusters(my_cdscomb, color_by = 'my_colourcom')

# Constructing GC lineage Single Cell Trajectories
expressed_genescomb <- row.names(subset(fData(my_cdscomb), num_cells_expressed >= 10))
my_cdscomb_subset <- my_cdscomb[expressed_genescomb, pData(my_cdscomb)$my_colourcom]
my_cdscomb_subset
my_cdscomb_subset <- detectGenes(my_cdscomb_subset, min_expr = 0.1)
fData(my_cdscomb_subset)$use_for_ordering <- fData(my_cdscomb_subset)$num_cells_expressed > 0.05 * ncol(my_cdscomb_subset)
table(fData(my_cdscomb_subset)$use_for_ordering)
plot_pc_variance_explained(my_cdscomb_subset, return_all = FALSE)
my_cdscomb_subset <- reduceDimension(my_cdscomb_subset, max_components = 2, num_dim = 20,
                                     reduction_method = 'tSNE', verbose = TRUE)
my_cdscomb_subset <- clusterCells(my_cdscomb_subset, num_clusters = 10)
table(pData(my_cdscomb_subset)$Cluster)

# Plots - GC lineage
plot_cell_clusters(my_cdscomb_subset)
plot_cell_clusters(my_cdscomb_subset, markers="Lgr5")
plot_cell_clusters(my_cdscomb_subset, markers="Esr1")
plot_cell_clusters(my_cdscomb_subset, markers="Esr2")

# DEG
clustering_DEG_genescomb <- differentialGeneTest(my_cdscomb_subset,
                                                 fullModelFormulaStr = '~Cluster',
                                                 cores = 4)
dim(clustering_DEG_genescomb)
clustering_DEG_genescomb %>% arrange(qval) %>% head()
my_ordering_genescomb <- row.names(clustering_DEG_genescomb)[order(clustering_DEG_genescomb$qval)][1:1000]
my_cdscomb_subset <- setOrderingFilter(my_cdscomb_subset, ordering_genes = my_ordering_genescomb)
my_cdscomb_subset <- reduceDimension(my_cdscomb_subset, method = 'DDRTree')
my_cdscomb_subset <- orderCells(my_cdscomb_subset)

## pseudotime

# differential GeneTest
head(pData(my_cdscomb_subset))
my_pseudotime_decomb <- differentialGeneTest(my_cdscomb_subset,
                                             fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                             cores = 8)
my_pseudotime_decomb %>% arrange(qval) %>% head()
my_pseudotime_decomb %>% arrange(qval) %>% head() %>% select(id) -> my_pseudotime_genecomb
my_pseudotime_genecomb <- my_pseudotime_genecomb$id
plot_genes_in_pseudotime(my_cdscomb_subset[my_pseudotime_genecomb,])

# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_decomb %>% arrange(qval) %>% head(50) %>% select(id) -> gene_to_clustercomb
gene_to_clustercomb <- gene_to_clustercomb$id
my_pseudotime_clustercomb <- plot_pseudotime_heatmap(my_cdscomb_subset[gene_to_clustercomb,],
                                                     num_clusters = 9,
                                                     cores = 8,
                                                     show_rownames = TRUE,
                                                     return_heatmap = TRUE)
newdatacomb <- data.frame(Pseudotime = seq(min(pData(my_cdscomb_subset)$Pseudotime), 
                                           max(pData(my_cdscomb_subset)$Pseudotime), length.out = 100))
my_clustercomb <- cutree(my_pseudotime_clustercomb$tree_row, 9)
my_clustercomb 

# genes in cluster 1
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 1]),"gene_short_name"]
# genes in cluster 2
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 2]),"gene_short_name"]
# genes in cluster 3
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 3]),"gene_short_name"]
# genes in cluster 4
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 4]),"gene_short_name"]
# genes in cluster 5
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 5]),"gene_short_name"]
# genes in cluster 6
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 6]),"gene_short_name"]
# genes in cluster 7
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 7]),"gene_short_name"]
# genes in cluster 8
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 8]),"gene_short_name"]
# genes in cluster 9
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 9]),"gene_short_name"]
# genes in cluster 10
my_pseudotime_decomb[names(my_clustercomb[my_clustercomb == 10]),"gene_short_name"]

# Set the root state (Lgr5+)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "State", markers="Lgr5", cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
my_cdscomb_subset <- orderCells(my_cdscomb_subset, root_state = 2)

# Plots
plot_cell_clusters(my_cdscomb_subset, color_by="Cluster")
plot_cell_clusters(my_cdscomb_subset, color_by="Pseudotime")
plot_cell_clusters(my_cdscomb_subset, color_by="State")
plot_cell_clusters(my_cdscomb_subset, markers="Lgr5")
plot_cell_clusters(my_cdscomb_subset, markers="Krt8")
plot_cell_clusters(my_cdscomb_subset, markers="Esr1")
plot_cell_clusters(my_cdscomb_subset, markers="Esr2")
plot_cell_clusters(my_cdscomb_subset, markers="Mki67")
plot_cell_clusters(my_cdscomb_subset, markers="Foxl2")
plot_cell_clusters(my_cdscomb_subset, markers="Col1a2")
plot_cell_clusters(my_cdscomb_subset, markers="Ramp2")
plot_cell_clusters(my_cdscomb_subset, markers="Fst")
plot_cell_clusters(my_cdscomb_subset, markers="Amh")
plot_cell_clusters(my_cdscomb_subset, markers="Inha")
plot_cell_clusters(my_cdscomb_subset, markers="Nr5a2")
plot_cell_clusters(my_cdscomb_subset, markers="Ihh")
plot_cell_clusters(my_cdscomb_subset, markers="Dhh")
plot_cell_clusters(my_cdscomb_subset, markers="Rspo1")
plot_cell_clusters(my_cdscomb_subset, markers="Tgfbr1")
plot_cell_clusters(my_cdscomb_subset, markers="Tgfbr2")
plot_cell_clusters(my_cdscomb_subset, markers="Tgfbr3")
plot_cell_clusters(my_cdscomb_subset, markers="Tgfb1")
plot_cell_clusters(my_cdscomb_subset, markers="Smad3")

plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Fst", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Wnt6", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Aldh1a1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Icam2", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Kctd14", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Lgals7", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Ecscr", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Emcn", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Igfbp7", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Gng11", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="S100acomb", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Cd34", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Esr1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, reverse = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Lgr5", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, reverse = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Foxl2", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, reverse = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Amh", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Inha", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Sycp3", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Krt8", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Lgals7", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Kazald1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Bex4", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Col27a1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Trim62", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Serpine2", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Vcan", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Inhbb", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Pax8", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Oct4", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
# Plots (PND7&14 Vs E16)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Dppa3", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Ooep", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Uchl1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Xdh", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Zp3", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Gtsf1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Nobox", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Oas1c", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Pinlyp", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
# Plots (E16 Vs PND7 and 14)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="H19", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Gng13", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Lgals7", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Smpx", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Ifitm1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Igf2", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Cst8", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Mest", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Podxl", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Rspo1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Tgfbr1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Tgfbr2", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Tgfbr3", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Tgfb1", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Mki67", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Smad3", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", markers="Fst", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 0)

plot_genes_branched_pseudotime(my_cdscomb_subset[my_gene,], branch_states = NULL, branch_point = 1,
                               branch_labels = NULL, method = "fitting", min_expr = NULL,
                               cell_size = 0.75, nrow = NULL, ncol = 1, panel_order = NULL,
                               color_by = "State", expression_curve_linetype_by = "Branch",
                               trend_formula = "~ sm.ns(Pseudotime, df=3) * Branch",
                               reducedModelFormulaStr = NULL, label_by_short_name = TRUE,
                               relative_expr = TRUE)

# cluster cell type
Krt18_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Krt18"))
Esr1_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Esr1"))
Esr2_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Esr2"))
Lgr5_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Lgr5"))
Foxl2_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Foxl2"))
Fst_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Fst"))
Ihh_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Ihh"))
Inha_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Inha"))
Wt1_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Wt1"))
Smad3_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Smad3"))
Amh_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Amh"))
Loxl1_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Loxl1"))
Col1a2_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Col1a2"))
Lgals7_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Lgals7"))
Ogn_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Ogn"))
Dcn_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Dcn"))
Ramp2_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Ramp2"))
Cyp11a1_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Cyp11a1"))
Cyp17a1_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Cyp17a1"))
Cyp19a1_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Cyp19a1"))
Plk2_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Plk2"))
Gli1_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Gli1"))
Mgp_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Mgp"))
Nr5a2_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Nr5a2"))
Rspo1_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Rspo1"))
Lgals7_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Lgals7"))
Nrg4_id <- row.names(subset(fData(my_cdscomb), gene_short_name == "Nrg4"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth07, "pre-GC", classify_func = function(x) { x[Foxl2_id,] >= 0.5 })
cth <- addCellType(cth07, "GCs(Ihh-)", classify_func = function(x) { x[Amh_id,] >= 0.5 }, parent_cell_type_name = "pre-GC")
cth <- addCellType(cth07, "GCs(Ihh+)", classify_func = function(x) { x[Ihh_id,] >= 0.5 }, parent_cell_type_name = "GC(Ihh-)")
cth <- addCellType(cth07, "OSE-Lgr5-", classify_func = function(x) { x[Krt8_id,] >= 0.5 & x[Foxl2_id,] < 0.5 })
cth <- addCellType(cth07, "OSE-Lgr5+", classify_func = function(x) { x[Lgr5_id,] >= 0.5 }, parent_cell_type_name = "OSE-Lgr5-")
my_cdscomb_subset <- classifyCells(my_cdscomb_subset, cth, 0.1)
table(pData(my_cdscomb_subset)$CellType)
pie <- ggplot(pData(my_cdscomb_subset), aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
disp_tablecomb <- dispersionTable(my_cdscomb_subset)
unsup_clustering_genescomb <- subset(disp_tablecomb, mean_expression >= 0.1)
my_cdscomb_subset <- setOrderingFilter(my_cdscomb_subset, unsup_clustering_genescomb$gene_id)
plot_ordering_genes(my_cdscomb_subset)
plot_pc_variance_explained(my_cdscomb_subset, return_all = F)
my_cdscomb_subset <- reduceDimension(my_cdscomb_subset, max_components = 2, num_dim = 10, reduction_method = 'tSNE', verbose = T)
my_cdscomb_subset <- clusterCells(my_cdscomb_subset, num_clusters = 6)

# Plots - cell type
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "Cluster")
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType")
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "Pseudotime")
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Esr1", "Krt18", "Lgr5", "Amh", "Esr2", "Inha"))
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Esr1", "Esr2", "Tgfbr1", "Tgfb1i1", "Lgr5", "Foxl2"))
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Esr1", "Esr2", "Tgfbr1", "Kitl", "Lgr5", "Foxl2"))
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Rspo1", "Fst", "Inha", "Krt8", "Cdkn1b", "Mki67"))
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Esr1", "Lgr5", "Krt8", "Fst", "Smad3", "Tgfbr1"))
my_cdscomb_subset <- reduceDimension(my_cdscomb_subset, max_components = 2, num_dim = 10, 
                                     reduction_method = 'tSNE', residualModelFormulaStr = "~num_genes_expressed", verbose = T)
my_cdscomb_subset <- clusterCells(my_cdscomb_subset, num_clusters = 6)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", cell_size = 3)
my_cdscomb_subset <- clusterCells(my_cdscomb_subset, num_clusters = 6)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "Cluster") + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Lgr5")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Amh")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Ihh")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Fst")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Esr2")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Esr1")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Foxl2")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Mki67")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Tgfb1")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Tgfbr1")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Tgfbr2")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Tgfbr3")) + facet_wrap(~CellType)
plot_cell_clusters(my_cdscomb_subset, 1, 2, color = "CellType", markers = c("Inha")) + facet_wrap(~CellType)

# DEG - Cell type
my_cdscomb_subset <- clusterCells(my_cdscomb_subset)
disp_tablecomb <- dispersionTable(my_cdscomb_subset)
ordering_genescomb <- subset(disp_tablecomb, mean_expression >= 0.1)
my_cdscomb_subset <- setOrderingFilter(my_cdscomb_subset, ordering_genescomb)
my_cdscomb_subset <- reduceDimension(my_cdscomb_subset)
my_cdscomb_subset<- orderCells(my_cdscomb_subset)
diff_test_rescomb <- differentialGeneTest(my_cdscomb_subset[expressed_genescomb,], fullModelFormulaStr = "~CellType")
ordering_genescomb <- row.names (subset(diff_test_rescomb, qval < 0.01))
my_cdscomb_order <- setOrderingFilter(my_cdscomb_subset, ordering_genescomb)
plot_ordering_genes(my_cdscomb_order)
my_cdscomb_order <- reduceDimension(my_cdscomb_order, max_components = 2, method = 'DDRTree')
my_cdscomb_order <- orderCells(my_cdscomb_order)

# Set the root state (Lgr5+)
my_cdscomb_order <- orderCells(my_cdscomb_order, root_state = 3)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", 
                     cell_size = 2, use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "State", 
                     cell_size = 2, use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "Pseudotime", 
                     cell_size = 2, use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", markers="Esr1", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", markers="Esr2", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", markers="Fst", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", markers="Lgr5", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", markers="Amh", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", markers="Inha", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", markers="Foxl2", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", markers="Krt18", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", markers="Mki67", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)

# plot_genes_in_pseudotime
GM_statecomb <- function(my_cdscomb_subset){
  if (length(unique(pData(my_cdscomb_subset)$State)) > 1){
    T0_counts <- table(pData(my_cdscomb_subset)$State, pData(my_cdscomb_subset)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "Pseudotime", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "State", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90) +
  facet_wrap(~State, nrow = 1)
plot_cell_trajectory(my_cdscomb_order, x = 1, y = 2, color_by = "CellType", 
                     use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90) +
  facet_wrap(~State, nrow = 1)
blast_genes <- row.names(subset(fData(my_cdscomb_order),
                                gene_short_name %in% c("Tgfbr1", "Lgr5", "Krt8", "Foxl2", "Esr2", "Nr5a1", "Esr1")))
plot_genes_jitter(my_cdscomb_order[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)
plot_genes_jitter(my_cdscomb_order[blast_genes,],
                  grouping = "CellType",
                  min_expr = 0.1)
combined_expressed_genes <-  row.names(subset(fData(my_cdscomb_order),
                                              num_cells_expressed >= 10))
combined_filtered <- my_cdscomb_order[combined_expressed_genes,]
combined_my_genes <- row.names(subset(fData(combined_filtered),
                                      gene_short_name %in% c("Lgr5", "Foxl2", "Esr2", "Esr1")))
combined_cds_subset <- combined_filtered[combined_my_genes,]
plot_genes_in_pseudotime(combined_cds_subset, color_by = "CellType", cell_size = 1.5)

branchTest(my_cdscomb_subset, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
           reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
           branch_states = NULL, branch_point = 1, relative_expr = TRUE,
           cores = 4, branch_labels = NULL, verbose = FALSE)

marker_diff <- markerDiffTable(my_cdscomb_subset[combined_expressed_genes,], cth, cores = 1)
semisup_clustering_genes <- row.names(marker_diff)[order(marker_diff$qval)][1:1000]
my_cdscomb_subset <- setOrderingFilter(my_cdscomb_subset, semisup_clustering_genes)
my_cdscomb_subset <- reduceDimension(my_cdscomb_subset, max_components = 2, method = 'DDRTree', norm_method = 'log')
my_cdscomb_subset <- orderCells(my_cdscomb_subset)
my_cdscomb_subset <- orderCells(my_cdscomb_subset, root_state = 3)
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "State", 
                     show_branch_points = FALSE, theta = 45) + theme(legend.position = "top")
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "CellType", 
                     show_branch_points = FALSE, theta = 45) + theme(legend.position = "top")
plot_cell_trajectory(my_cdscomb_subset, x = 1, y = 2, color_by = "Pseudotime", 
                     show_branch_points = FALSE, theta = 45) + theme(legend.position = "top")
combined_filtered <- my_cdscomb_subset[combined_expressed_genes,]
my_genes <- row.names(subset(fData(combined_filtered), gene_short_name %in% c("Lgr5", "Foxl2", "Esr2", "Esr1")))

cds_subset <- combined_filtered[my_genes,]

plot_genes_branched_pseudotime(cds_subset, branch_point = 1, color_by = "CellType", ncol = 1, cell_size = 1.5)
plot_genes_branched_pseudotime(cds_subset, branch_point = 1, color_by = "Pseudotime", ncol = 1, cell_size = 1.5)
marker_genes <- row.names(subset(fData(my_cdscomb_subset), gene_short_name %in% c("Lgr5", "Foxl2", "Esr2", "Esr1")))
diff_test_res <- differentialGeneTest(my_cdscomb_subset[marker_genes,], fullModelFormulaStr = "~CellType")

# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes[,c("gene_short_name", "pval", "qval")]
diff_test_res <- differentialGeneTest(my_cdscomb_order[marker_genes,], fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(my_cdscomb_order[sig_gene_names,],
                        num_clusters = 8,
                        cores = 4,
                        show_rownames = T)

# differential GeneTest
to_be_tested <- row.names(subset(fData(my_cdscomb_subset), 
                                 gene_short_name %in% c("Esr1", "Esr2", "Foxl2", "Lgr5", "Fst", "Wt1", "Gli1")))
cds_subset <- my_cdscomb_subset[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~CellType")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_jitter(cds_subset,
                  grouping = "CellType",
                  color_by = "CellType",
                  nrow= 1,
                  ncol = NULL,
                  plot_trend = TRUE)
full_model_fits <- fitModel(cds_subset,  modelFormulaStr = "~CellType")
reduced_model_fits <- fitModel(cds_subset, modelFormulaStr = "~1")
diff_test_res <- compareModels(full_model_fits, reduced_model_fits)
diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_in_pseudotime(cds_subset, color_by = "CellType", cell_size = 1.5)
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime", cell_size = 1.5)




