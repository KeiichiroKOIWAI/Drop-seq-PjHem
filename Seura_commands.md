# read libraries
```{r}
library("dplyr")
library("Seurat")
library("SeuratData")
library("Matrix")
library("tidyverse")
library("patchwork")
library("readr")
```
# read STARsolo filtered result
```{r}
X.mtx <- readMM("matrix.mtx")
rownames(X.mtx) <- read_tsv("features.tsv", col_names=FALSE)[, 2, drop=TRUE]
colnames(X.mtx) <- read_tsv("barcodes.tsv", col_names=FALSE)[, 1, drop=TRUE]

X <- CreateSeuratObject(counts=X.mtx, project = "X", names.field = 1)

# Create a new metadata column called "original" with the contents of "orig.ident"
X$original <- X$orig.ident
# Remove the "orig.ident" metadata column
X$orig.ident <- NULL
X@meta.data$"orig.ident" <- X@project.name
X <- SetIdent(X, value = "orig.ident")
X$stim <- "cont" or "infected"
X$sex <- "female" or "male"
X$protocol <- "fix&4+11" or "unfixed&4+12"
X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = "^MT-")
X[["percent.WSSV"]] <- PercentageFeatureSet(X, pattern = "^SWSSV-")
```
```{r}
# Visualize QC metrics as a violin plot
VlnPlot(X, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.WSSV"), ncol = 2)
```
# extract single cell data
```{r}
Xe <- subset(X, subset = nCount_RNA > 500 & nCount_RNA < 4000 &  percent.mt < 10)
VlnPlot(Xe, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.WSSV"), ncol = 2)
```

=================================================================================================
# SCTransform
```{r}
Xe <- SCTransform(Xe, vst.flavor = "v2", verbose = FALSE, vars.to.regress = c("percent.mt","percent.WSSV")) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.5, verbose = FALSE)
```
```{r}
p1 <- DimPlot(Xe, label = T, repel = T) + ggtitle("Unsupervised clustering")
p1
```
# Perform integration using pearson residuals
```{r}
Ye <- SCTransform(Ye, vst.flavor = "v2", verbose = FALSE, vars.to.regress = c("percent.mt","percent.WSSV")) %>%
    RunPCA(npcs = 30, verbose = FALSE)
```

# selecting a list of informative features
```{r}
list = list(Xe, Ye, Ze....)
# Choose the features to use when integrating multiple datasets. This function ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. It returns the top scoring features by this ranking.
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
```
# FindIntegrationAnchors
```{r}
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
```
# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells:
```{r}
combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
```
```{r}
combined.sct <- FindClusters(combined.sct, resolution = 0.5)
# To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.

DimPlot(combined.sct, reduction = "umap", split.by = "stim")
```
```{r}
p1 <- DimPlot(combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
p1 | p2
```
# add metadata "HemX"
```{r}
Idents(combined.sct) <- "seurat_clusters"
new.cluster.ids <- c("Hem4", "Hem6", "Hem2", "Hem5", "Hem1", "Hem7", "Hem3", "Hem8")
names(new.cluster.ids) <- levels(combined.sct)
combined.sct <- RenameIdents(combined.sct, new.cluster.ids)
combined.sct@active.ident <- factor(combined.sct@active.ident, levels=c("Hem1", "Hem2", "Hem3", "Hem4", "Hem5", "Hem6", "Hem7", "Hem8"))
combined.sct$hem_clusters <- combined.sct@active.ident
```
```{r}
combined.sct$celltype.stim <- paste(combined.sct$hem_clusters, combined.sct$stim, sep = "_")
```

# Identify differential expressed genes across conditions
# Using the normalized datasets with known celltype annotation, we can ask what genes change in different conditions for cells of the same type. First, we create a column in the meta.data slot to hold both the cell type and stimulation information and switch the current ident to that column.
# add meta data
# Find markers
```{r}
# To run differential expression, we make use of ‘corrected counts’ that are stored in the data slot of the the SCT assay. Corrected counts are obtained by setting the sequencing depth for all the cells to a fixed value and reversing the learned regularized negative-binomial regression model. Prior to performing differential expression, we first run PrepSCTFindMarkers, which ensures that the fixed value is set properly. Then we use FindMarkers(assay="SCT") to find differentially expressed genes. Here, we aim to identify genes that are differently expressed between stimulated and control B cells.
DefaultAssay(combined.sct) <- "integrated"
combined.sct <- PrepSCTFindMarkers(combined.sct)
Idents(combined.sct) <- "stim"
wssv.response <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "cont", ident.2 = "wssv", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
Idents(combined.sct) <- "celltype.stim"
wssv.response.Hem1 <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "Hem1_cont", ident.2 = "Hem1_wssv", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
wssv.response.Hem2 <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "Hem2_cont", ident.2 = "Hem2_wssv", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
wssv.response.Hem3 <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "Hem3_cont", ident.2 = "Hem3_wssv", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
wssv.response.Hem4 <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "Hem4_cont", ident.2 = "Hem4_wssv", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
wssv.response.Hem5 <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "Hem5_cont", ident.2 = "Hem5_wssv", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
wssv.response.Hem6 <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "Hem6_cont", ident.2 = "Hem6_wssv", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
wssv.response.Hem7 <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "Hem7_cont", ident.2 = "Hem7_wssv", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
wssv.response.Hem8 <- FindMarkers(combined.sct, assay = "SCT", ident.1 = "Hem8_cont", ident.2 = "Hem8_wssv", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
DefaultAssay(combined.sct) <- "SCT"
```
# Identify conserved cell type markers
# To identify canonical cell type marker genes that are conserved across conditions, we provide the FindConservedMarkers() function. This function performs differential gene expression testing for each dataset/group and combines the p-values using meta-analysis methods from the MetaDE R package. For example, we can identify genes that are conserved markers irrespective of stimulation condition in NK cells. Note that the PrepSCTFindMarkers command does not to be rerun here.
```{r}
DefaultAssay(combined.sct) <- "integrated"
Idents(combined.sct) <- "hem_clusters"
Hem1.markers <- FindConservedMarkers(combined.sct, assay = "SCT", ident.1 = "Hem1", grouping.var = "stim", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
Hem2.markers <- FindConservedMarkers(combined.sct, assay = "SCT", ident.1 = "Hem2", grouping.var = "stim", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
Hem3.markers <- FindConservedMarkers(combined.sct, assay = "SCT", ident.1 = "Hem3", grouping.var = "stim", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
Hem4.markers <- FindConservedMarkers(combined.sct, assay = "SCT", ident.1 = "Hem4", grouping.var = "stim", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
Hem5.markers <- FindConservedMarkers(combined.sct, assay = "SCT", ident.1 = "Hem5", grouping.var = "stim", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
Hem6.markers <- FindConservedMarkers(combined.sct, assay = "SCT", ident.1 = "Hem6", grouping.var = "stim", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
Hem7.markers <- FindConservedMarkers(combined.sct, assay = "SCT", ident.1 = "Hem7", grouping.var = "stim", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
Hem8.markers <- FindConservedMarkers(combined.sct, assay = "SCT", ident.1 = "Hem8", grouping.var = "stim", min.pct = 0.5, logfc.threshold = 0.5, verbose = FALSE)
DefaultAssay(combined.sct) <- "SCT"
```
