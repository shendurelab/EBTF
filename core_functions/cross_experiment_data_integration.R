
######################################
### Cross-experiment data integration 
### contact: cxqiu@uw.edu

#########################################################
### Step-1: Integration of multiple datasets using Seurat

print("Loading Monocle3 and Seurat")
suppressMessages(library(monocle3))
suppressMessages(library(Seurat))

print("Loading packages for regular data analysis, e.g. dplyr")
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))

print("Loading packages for plotting, e.g. ggplot2")
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(gridExtra))
suppressMessages(library(viridis))
suppressMessages(library(gplots))

library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 80000 * 1024^2)   ### adjust the parameter based on the size of datasets

### Each dataset in Seurat object

obj_1 = readRDS("obj_1.rds")
count_1 = GetAssayData(obj_1, slot = "counts")
pd_1 = data.frame(obj_1[[]])
pd_1$group = "dataset_1"

obj_2 = readRDS("obj_2.rds")
count_2 = GetAssayData(obj_2, slot = "counts")
pd_2 = data.frame(obj_2[[]])
pd_2$group = "dataset_2"

gene_overlap = intersect(rownames(count_1), rownames(count_2))
col_include = intersect(colnames(pd_1), colnames(pd_2))

obj_1 = CreateSeuratObject(count_1[gene_overlap,], meta.data = pd_1[,col_include])
obj_2 = CreateSeuratObject(count_2[gene_overlap,], meta.data = pd_2[,col_include])

obj = merge(obj_1, obj_2)

print(table(obj$group))
print(dim(obj))

obj.list <- SplitObject(obj, split.by = "group")
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = "rpca", 
                                  dims = 1:50)
obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, min.dist = 0.3, n.components = 2)

saveRDS(obj.integrated, "obj_integration.rds")

obj.integrated$UMAP_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
obj.integrated$UMAP_2 = Embeddings(obj.integrated, reduction = "umap")[,2]
saveRDS(data.frame(obj.integrated[[]]), "obj_integration_pd.rds")


##################################################
### Step-2: Visualizing integration result by UMAP

p1 = ggplot(pd) +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
    geom_point(data = subset(pd, group == "dataset_1"),
               aes(x = UMAP_1, y = UMAP_2), color = "red", size=0.15) +
    theme_void() +
    theme(legend.position="none")

p2 = ggplot(pd) +
    geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
    geom_point(data = subset(pd, group == "dataset_2"),
               aes(x = UMAP_1, y = UMAP_2), color = "red", size=0.15) +
    theme_void() +
    theme(legend.position="none")

ggsave("UMAP_integration.png", p1 + p2, dpi = 300, height = 5, width = 10)



