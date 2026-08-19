# loading packages
library(Seurat)
library(Matrix)

# Importing data
data_dir <- "C:/Users/ommiv/Downloads/GSE134520_RAW(1)"
# use forward slash here
list.files(data_dir)

# defining the samples
sample_files <- c(
  IMS1 = "GSM3954954_processed_IMS1.txt.gz",
  IMS2 = "GSM3954955_processed_IMS2.txt.gz",
  IMS3 = "GSM3954956_processed_IMS3.txt.gz",
  IMS4 = "GSM3954957_processed_IMS4.txt.gz",
  IMW1 = "GSM3954952_processed_IMW1.txt.gz",
  IMW2 = "GSM3954953_processed_IMW2.txt.gz",
  CAG1 = "GSM3954949_processed_CAG1.txt.gz",
  CAG2 = "GSM3954950_processed_CAG2.txt.gz",
  CAG3 = "GSM3954951_processed_CAG3.txt.gz",
  NAG1 = "GSM3954946_processed_NAG1.txt.gz",
  NAG2 = "GSM3954947_processed_NAG2.txt.gz",
  NAG3 = "GSM3954948_processed_NAG3.txt.gz",
  EGC = "GSM3954958_processed_EGC.txt.gz"
)

sample_files

file.exists(test_file)
test_file <- file.choose()
test_file
file.exists(test_file)

# reading the IMS1 file
test_data <- read.delim(
  test_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
dim(test_data)
# gives you genes, cells
head(test_data[, 1:5])
head(rownames(test_data)) # genes
head(colnames(test_data)) # cells/barcodes

# converting data to a sparse matrix
# until now we had a R data frame(test_data)
# seurat works with sparse matrices bcz most entries in Sc matrix are 0

# checking data type
class(test_data)

# converting to matrix
test_matrix <- as.matrix(test_data)
class(test_matrix)

# converting to sparse matrix
library(Matrix)
test_matrix <- Matrix(
  test_matrix,
  sparse = TRUE
)
class(test_matrix)
# check dimensions of matrix
dim(test_matrix)
head(rownames(test_matrix))
head(colnames(test_matrix))

# checking for duplicated genes
sum(duplicated(rownames(test_matrix)))
# if ans = 0, no duplicated gene names

# checking for missing values
sum(is.na(test_matrix))

# checking for count values
summary(as.numeric(test_matrix))
# this will give min, max, median...
# checking values are integers
all(as.numeric(test_matrix)) == floor(as.numeric(test_matrix))

class(test_matrix)
dim(test_matrix)
sum(duplicated(rownames(test_matrix)))
sum(is.na(test_matrix))

# ___________

# creating seurat package 
library(Seurat)
IMS1 <- CreateSeuratObject(
  counts = test_matrix,
  project = "GSE134520",
  min.cells = 3,
  min.features = 200,
) 

IMS1
head(IMS1@meta.data)

# QC

IMS1[["percent.mt"]] <- PercentageFeatureSet(IMS1, pattern = "^MT-")

# Qc metrics as violin plot
VlnPlot(IMS1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)

# FeatureScatter 
plot1 <- FeatureScatter(IMS1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(IMS1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

IMS1 <- subset(IMS1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <5)

#___________

# Normalizing the data
IMS1 <- NormalizeData(IMS1, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features(feature selection)
IMS1 <- FindVariableFeatures(IMS1, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(IMS1), 10)

# plotting
plot1 <- VariableFeaturePlot(IMS1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
all.genes <- rownames(IMS1)
IMS1 <- ScaleData(IMS1, features = all.genes)

# Performing linear dimensional reduction
IMS1 <- RunPCA(IMS1, features = VariableFeatures(object = IMS1))

print(IMS1[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(IMS1, dims = 1:2, reduction ="pca")

DimPlot(IMS1, reduction = "pca") + NoLegend()

DimHeatmap(IMS1, dims =1, cells = 500, balanced = TRUE)

DimHeatmap(IMS1, dims = 1:15, cells = 500, balanced = TRUE)

# determine the dimensionality of dataset
ElbowPlot(IMS1)

# clustering the cells
IMS1 <- FindNeighbors(IMS1, dims = 1:10)
IMS1 <- FindClusters(IMS1, resolution = 0.5)

head(Idents(IMS1), 5) # cluster Id of first 5 cells

# Non linear dimensionality reduction
IMS1 <- RunUMAP(IMS1,dims = 1:10)

DimPlot(IMS1, reduction = "umap")

saveRDS(IMS1, file = "C:/Users/ommiv/Downloads/IMS1_seurat.rds")

# finding cluster biomarkers
# finding all markers of cluster 2
cluster2.markers <- FindMarkers(IMS1, ident.1 = 2)
head(cluster2.markers, n = 5)

# finding all markers distinguishing cluster 5 from cluster 0 & 3
cluster5.markers <- FindMarkers(IMS1, ident.1 = 5, ident.2 = c(0,3))
head(cluster5.markers, n = 5)

# compared to all remaining cells
IMS1.markers <- FindAllMarkers(IMS1, only.pos = TRUE)
IMS1.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >1)

cluster0.markers <- FindMarkers(IMS1, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(IMS1, features = c("RAMP3", "PLVAP"))

# you can plot raw counts as well
VlnPlot(IMS1, features = c("KDR", "FABP5"), slot = "counts", log = TRUE)

FeaturePlot(IMS1, features = c("RAMP3", "ESAM", "CD320", "CD93", "VWF", "PODXL", "CD34", "VWA1",
                               "IPO11"))

IMS1.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(IMS1, features = top10$gene) + NoLegend()

