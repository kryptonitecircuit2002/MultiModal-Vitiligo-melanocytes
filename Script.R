library(dplyr)
library(Seurat)
library(patchwork)

#Create lists 
datapaths <- list(
  ".../PtA_processed",
  "..../PtC_processed",
  ".../PtD_processed",
  ".../PtE_processed",
  ".../PtF_processed",
  ".../PtG_processed"
)
# Create a list to store Seurat objects
data_list <- lapply(datapaths, function(path) {
  data <- Read10X(data.dir = path)
  combined_patient_obj <- CreateSeuratObject(counts = data, project = basename(path))
  return(combined_patient_obj)
})

all_patient_data <- merge(data_list[[1]], y = data_list[-1])

# Strip off the patient-specific suffix (_1, _2, etc.) 
original_barcodes <- gsub("_\\d+$", "", colnames(all_patient_data))

# Create metadata for seurat indicating whether the barcode belongs to lesional or non-lesional tissue
all_patient_data$sampletype  <- ifelse(grepl("-1$", original_barcodes), "non-lesional", "lesional")

#split the data by sample type
all_patient_data[["RNA"]] <- split(all_patient_data[["RNA"]], f = all_patient_data$sampletype)

#Perform QC
all_patient_data[["percent.mt"]] <- PercentageFeatureSet(all_patient_data, pattern = "^MT-")
VlnPlot(all_patient_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Filter samples
all_patient_data <- subset(all_patient_data, subset = nFeature_RNA < 7000 & nCount_RNA < 40000 & percent.mt < 10)

#Calculate number of cells across both 
table(all_patient_data$sampletype)


-------------------------------------------------------------------Downsampling-----------------------------------------------------------------------------
  
#Downsample them 
# Get all cell names
all_cells <- colnames(all_patient_data)

# Sample lesional cells
lesional_cells <- all_cells[all_patient_data$sampletype == "lesional"]
non_lesional_cells <- all_cells[all_patient_data$sampletype == "non-lesional"]

# Sample 7800 cells from each group, if available
set.seed(123)  # For reproducibility
sampled_lesional <- if (length(lesional_cells) >= 7800) {
  sample(lesional_cells, 7800)
} else {
  lesional_cells
}

sampled_non_lesional <- if (length(non_lesional_cells) >= 7800) {
  sample(non_lesional_cells, 7800)
} else {
  non_lesional_cells
}

# Combine sampled cells
final_sampled_cells <- c(sampled_lesional, sampled_non_lesional)

# Subset the original Seurat object
all_patient_data.subsampled <- subset(all_patient_data, cells = final_sampled_cells)

#Calculate number of cells across both 
table(all_patient_data.subsampled$sampletype)


# run standard anlaysis workflow
all_patient_data.subsampled <- NormalizeData(all_patient_data.subsampled)
all_patient_data.subsampled <- FindVariableFeatures(all_patient_data.subsampled)
all_patient_data.subsampled <- ScaleData(all_patient_data.subsampled)
all_patient_data.subsampled <- RunPCA(all_patient_data.subsampled)

all_patient_data.subsampled <- FindNeighbors(all_patient_data.subsampled, dims = 1:30, reduction = "pca")
all_patient_data.subsampled <- FindClusters(all_patient_data.subsampled, resolution = 1, cluster.name = "unintegrated_clusters")
all_patient_data.subsampled <- RunUMAP(all_patient_data.subsampled, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(all_patient_data.subsampled, reduction = "umap.unintegrated", group.by = c("sampletype", "seurat_clusters"))
DimPlot(all_patient_data.subsampled, reduction = "umap.unintegrated", split.by = "sampletype", label = TRUE)


-------------------------------------------------------------------Integration-------------------------------------------------------------------------------------------------------------
#Perform Integration
all_patient_data.integrated <- IntegrateLayers(object = all_patient_data.subsampled, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

# re-join layers after integration
all_patient_data.integrated[["RNA"]] <- JoinLayers(all_patient_data.integrated[["RNA"]])

all_patient_data.integrated <- FindNeighbors(all_patient_data.integrated, reduction = "integrated.cca", dims = 1:30)
all_patient_data.integrated <- FindClusters(all_patient_data.integrated, resolution = 0.8)
all_patient_data.integrated <- RunUMAP(all_patient_data.integrated, dims = 1:30, reduction = "integrated.cca")
DimPlot(all_patient_data.integrated, reduction = "umap", group.by = c("sampletype", "seurat_clusters"))
DimPlot(all_patient_data.integrated, reduction = "umap", split.by = "sampletype", label = TRUE)

#Check for Melanocyte cluster
DotPlot(all_patient_data.integrated, features = c("ALG3", "SUDS3", "ORAI1", "EIF4BP1", "S100B", "CFAP36", "GPS1", "CNP"), split.by = "sampletype")+ 
  scale_color_viridis(option = "G")
DoHeatmap(all_patient_data.integrated,
          features = c("ALG3", "SUDS3", "ORAI1", "EIF4BP1", "S100B", "CFAP36", "GPS1", "CNP"),
          group.by = "seurat_clusters",
          label = TRUE) + 
  scale_fill_viridis(option = "A")

#Find conserved markers
melanocyte.markers <- FindConservedMarkers(all_patient_data.integrated, ident.1 = 21, grouping.var = "sampletype", verbose = FALSE)
write.csv(melanocyte.markers, file = "cluster21_markers_results.csv", row.names = TRUE)

#Perform DEG
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

aggregate_data <- AggregateExpression(all_patient_data.integrated, group.by = c("seurat_clusters", "sampletype"), return.seurat = TRUE)
VlnPlot(all_patient_data.integrated, features = "CD45", split.by = "sampletype", group.by = "seurat_clusters",
        pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

top10 <- melanocyte.markers %>%
  mutate(sampletype = if_else(lesional_avg_log2FC > 'non-lesional_avg_log2FC', "lesional", "non-lesional")) %>%
  slice_max(order_by = if_else(sampletype == "lesional", lesional_avg_log2FC, non-lesional_avg_log2FC), n = 10)

RidgePlot(all_patient_data.integrated, features = c("ALG3", "SUDS3", "ORAI1", "EIF4BP1", "S100B", "CFAP36", "GPS1", "CNP"))

