
######################################################
### Merging samples and performing dimension reduction
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

run_list = paste0("E2E_rep", c("1_1", "1_2"))

count = NULL
pd = NULL
for(run_id in run_list){
    print(run_id)
    obj_tmp = readRDS(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/obj_monoclonal_EB_proof_", run_id, ".rds"))
    print(dim(obj_tmp))
    
    count = cbind(count, GetAssayData(obj_tmp, slot = "counts"))
    pd = rbind(pd, data.frame(obj_tmp[[]]))
}

saveRDS(count, paste0(work_path, "/analysis/monoclonal_EB_proof/backup/count.rds"))
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_proof/backup/pd.rds"))


######################################################################
### Step-1: performing regular seurat analysis but including align_cds

count = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/backup/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/backup/pd.rds"))

### only include rep_1_1 and rep_1_2
keep = pd$experiment_id %in% paste0("monoclonal_EB_proof_E2E_rep", c("1_1", "1_2"))
count = count[,keep]
pd = pd[keep,]

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]

obj = CreateSeuratObject(count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
obj = ScaleData(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 30, verbose = FALSE)

### after calculating PCA, performing align_cds to correct lanes
### pca_coor = as.matrix(Embeddings(obj, reduction = "pca"))
### aligned_coor = batchelor::reducedMNN(as.matrix(pca_coor),
###                                      batch = data.frame(obj[[]])[,"experiment_id"],
###                                      k=20)
### aligned_coor = aligned_coor$corrected
### 
### colnames(aligned_coor) = paste0("aligned_", 1:30)
### rownames(aligned_coor) = colnames(obj)
### 
### obj[['aligned']] = Seurat::CreateDimReducObject(embeddings=as.matrix(aligned_coor), key='aligned_')

### then, we perform UMAP or clustering based on the corrected PCA coordinates
obj = RunUMAP(object = obj, reduction = "pca", dims = 1:30, min.dist = 0.3, verbose = FALSE)
obj = FindNeighbors(object = obj, reduction = "pca", dims = 1:30, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 1, verbose = FALSE)
saveRDS(obj, paste0(work_path, "/analysis/monoclonal_EB_proof/backup/obj_aligned_first_try.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_proof/backup/obj_aligned_pd_first_try.rds"))


###################################################################
### Step-2: removing some clusters with very low UMI and high RIBO%

count = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/backup/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/backup/pd.rds"))

### only include rep_1_1 and rep_1_2
keep = pd$experiment_id %in% paste0("monoclonal_EB_proof_E2E_rep", c("1_1", "1_2"))
pd = pd[keep,]
count = count[,keep]

### In the first try, I found there are several cell clusters showing lower UMI 
pd_x = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/backup/obj_aligned_pd_first_try.rds"))
print(sum(rownames(pd) == rownames(pd_x)))
keep = !pd_x$RNA_snn_res.1 %in% c(14,20)

pd = pd[keep,]
count = count[,keep]

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]

### regular processing
obj = CreateSeuratObject(count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
obj = ScaleData(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 30, verbose = FALSE)

### then, we perform UMAP or clustering based on the corrected PCA coordinates
obj = RunUMAP(object = obj, reduction = "pca", dims = 1:30, min.dist = 0.3, verbose = FALSE)
obj = FindNeighbors(object = obj, reduction = "pca", dims = 1:30, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 1, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 2, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 3, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 5, verbose = FALSE)
saveRDS(obj, paste0(work_path, "/analysis/monoclonal_EB_proof/obj_processed.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_proof/obj_processed_pd.rds"))


#################################
### Step-3: cell type annotations

anno = rep(NA, nrow(pd))
anno[pd$RNA_snn_res.1 %in% c(10)] = "Primordial germ cells"
anno[pd$RNA_snn_res.1 %in% c(4,6)] = "Epiblast"
anno[pd$RNA_snn_res.1 %in% c(11)] = "Primitive streak"
anno[pd$RNA_snn_res.1 %in% c(7,13)] = "Neuroectoderm"
anno[pd$RNA_snn_res.1 %in% c(14)] = "Early neurons"
anno[pd$RNA_snn_res.1 %in% c(8)] = "Surface ectoderm"
anno[pd$RNA_snn_res.1 %in% c(12)] = "Definitive endoderm"
anno[pd$RNA_snn_res.1 %in% c(16)] = "Extraembryonic visceral endoderm"
anno[pd$RNA_snn_res.1 %in% c(17)] = "Nascent mesoderm"
anno[pd$RNA_snn_res.1 %in% c(1,2)] = "Lateral plate mesoderm"
anno[pd$RNA_snn_res.1 %in% c(3)] = "Cardiomyocytes"
anno[pd$RNA_snn_res.1 %in% c(0,5)] = "Paraxial mesoderm" ### Meox1, Meox2, Tbx1, Pax3
anno[pd$RNA_snn_res.1 %in% c(9)] = "Endothelial cells"
anno[pd$RNA_snn_res.1 %in% c(15)] = "Erythroid cells"
anno[pd$RNA_snn_res.1 %in% c(18)] = "Blood progenitors"

anno[pd$RNA_snn_res.5 %in% c(45)] = "Notochord"
anno[pd$RNA_snn_res.5 %in% c(31)] = "Neuroectoderm"
anno[pd$RNA_snn_res.1 %in% c(1,2) & !is.na(pd_x$celltype) & pd_x$celltype == "Cardiomyocytes"] = "Cardiomyocytes"

pd$celltype = as.vector(anno)
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_proof/obj_processed_pd.rds"))
### n = 19,059 cells

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.3, color = "black") +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.15) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/monoclonal_EB_proof_celltype.png"),
               dpi = 300,
               height  = 4, 
               width = 5), silent = TRUE)


try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.3, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.15) +
        ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/monoclonal_EB_proof_celltype_label.png"),
               dpi = 300,
               height  = 4, 
               width = 5), silent = TRUE)






