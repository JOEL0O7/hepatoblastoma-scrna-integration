# script to integrate scRNA-Seq datasets to correct for batch effects
# data: Modeling Hepatoblastoma: Identification of Distinct Tumor Cell Populations and Key Genetic Mechanisms through Single Cell Sequencing (scRNA-seq)
# data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180665

setwd("/home/joel/Bioinformatics/Seurat practice/Integrate scRNA Dataset")

# Load Libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# get data location
dirs <- list.dirs( recursive = F, full.names = F)

# Create Seurat obj for each directories
for (x in dirs){
  name <- gsub("_filtered_feature_bc_matrix","",x)  # Clean the name
  
  # Read 10X matrix
  mtx_obj <- ReadMtx(mtx = paste0(x,"/matrix.mtx.gz"),
                     features = paste0(x,"/features.tsv.gz"),
                     cells = paste0(x,"/barcodes.tsv.gz"))
  
  # Creating a Seurat obj and assigning it to cleaned name      
  seurat_obj <- CreateSeuratObject(counts = mtx_obj)
  assign(name, seurat_obj)
}

# merge datasets
  
    # It would be tedious and time consuming to perform QC on each seurat obj one by one, so we merge all into one and then perform QC
ls()
merged_seurat <- merge(HB17_background, y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor, HB53_background,HB53_tumor),
      add.cell.ids = ls()[2:8], 
      project = "HB")
merged_seurat



  #mat <- GetAssayData(seurat_obj, slot= "counts")
  #mat[1:50,6:10]

# QC and filtering----------------
View(merged_seurat@meta.data)

# Create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# Split the sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = "sample", into = c("Patient", "Type", "Barcode"), sep = "_")

# Calculate mitochondrial percentage
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")

# Explore QC
hist(merged_seurat$nFeature_RNA, breaks = 50, main = "nFeature_RNA distribution")
hist(merged_seurat$nCount_RNA, breaks = 50, main = "nCount_RNA distribution")
hist(merged_seurat$percent.mt, breaks = 50, main = "npercent.mt")

FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")


# Filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 & nFeature_RNA > 500 & percent.mt < 10)
merged_seurat_filtered


# perform standard workflow steps to figure out if we see any batch effects --------

merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)                # Normalization
merged_seurat_filtered <- FindVariableFeatures(merged_seurat_filtered)         # Identify (default = 2000) highly variable features



merged_seurat_filtered <- ScaleData(merged_seurat_filtered)                    # Scaling
merged_seurat_filtered <- RunPCA(merged_seurat_filtered)                       # Linear Dimensionality Reduction - PCA
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20)         # Non-linear Dimensionality Reduction - UMAP
DimPlot(merged_seurat_filtered, reduction = 'umap', label = T)


# Plot
p1 <- DimPlot(merged_seurat_filtered, reduction = "umap", group.by = "Patient")
p2 <- DimPlot(merged_seurat_filtered, reduction = "umap", group.by = "Type", cols = c("orange",'purple','brown'))

grid.arrange(p1,p2, ncol = 2)


# Perform integration to correct for batch effect
obj.list <- SplitObject(merged_seurat_filtered, split.by = "Patient")

for (i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
  
  # Collapse multiple layers if present (v5 fix)
  # This resolves the 'GetAssayData doesn't work for multiple layers' error.
  obj.list[[i]] <- JoinLayers(obj.list[[i]]) 
}

                

# Select integration feature
features <- SelectIntegrationFeatures(object.list = obj.list)                          #Stores the common variable features across 3 patients 

# Find Integration Anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features ) #Cells that express these features.

# Integrate Data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale, run PCA and UMAP on integrated data
seurat.integrated <- ScaleData(seurat.integrated)
seurat.integrated <- RunPCA(seurat.integrated)
seurat.integrated <- RunUMAP(seurat.integrated, dims = 1:50)


p3 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Patient")
p4 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Type", cols = c("orange",'purple','brown'))
grid.arrange(p3,p4, ncol = 2)

grid.arrange(p1,p2,p3,p4, nrow = 2, ncol = 2)

