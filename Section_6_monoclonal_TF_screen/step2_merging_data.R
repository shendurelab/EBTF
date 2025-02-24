
######################################################
### Merging samples and performing dimension reduction
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

run_list = paste0("clonal_EBTF_", c(1:16))

count = NULL
pd = NULL
for(run_id in run_list){
    print(run_id)
    obj_tmp = readRDS(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/obj_monoclonal_EB_TF_screen_", run_id, ".rds"))
    print(dim(obj_tmp))
    
    count = cbind(count, GetAssayData(obj_tmp, slot = "counts"))
    pd = rbind(pd, data.frame(obj_tmp[[]]))
}

saveRDS(count, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/backup/count.rds"))
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/backup/pd.rds"))


######################################################################
### Step-1: performing regular seurat analysis but including align_cds

source("~/work/scripts/utils.R")
mouse_gene <- read.table("~/work/tome/code/mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

count = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/backup/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/backup/pd.rds"))

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]

obj = CreateSeuratObject(count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
obj = ScaleData(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 30, verbose = FALSE)

### after calculating PCA, performing align_cds to correct lanes
 pca_coor = as.matrix(Embeddings(obj, reduction = "pca"))
 aligned_coor = batchelor::reducedMNN(as.matrix(pca_coor),
                                      batch = data.frame(obj[[]])[,"experiment_id"],
                                      k=20)
 aligned_coor = aligned_coor$corrected
 
 colnames(aligned_coor) = paste0("aligned_", 1:30)
 rownames(aligned_coor) = colnames(obj)
 
 obj[['aligned']] = Seurat::CreateDimReducObject(embeddings=as.matrix(aligned_coor), key='aligned_')

### then, we perform UMAP or clustering based on the corrected PCA coordinates
obj = RunUMAP(object = obj, reduction = "aligned", dims = 1:30, min.dist = 0.3, verbose = FALSE)
obj = FindNeighbors(object = obj, reduction = "aligned", dims = 1:30, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 1, verbose = FALSE)
saveRDS(obj, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/backup/obj_aligned_first_try.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/backup/obj_aligned_pd_first_try.rds"))


###################################################################
### Step-2: removing some clusters with very low UMI and high RIBO%


count = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/backup/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/backup/pd.rds"))

### In the first try, I found there are several cell clusters showing lower UMI 
pd_x = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/backup/obj_aligned_pd_first_try.rds"))
keep = !pd_x$RNA_snn_res.1 %in% c(18)

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

### after calculating PCA, performing align_cds to correct lanes
  pca_coor = as.matrix(Embeddings(obj, reduction = "pca"))
  aligned_coor = batchelor::reducedMNN(as.matrix(pca_coor),
                                    batch = data.frame(obj[[]])[,"experiment_id"],
                                    k=20)
  aligned_coor = aligned_coor$corrected
  
  colnames(aligned_coor) = paste0("aligned_", 1:30)
  rownames(aligned_coor) = colnames(obj)
  
  obj[['aligned']] = Seurat::CreateDimReducObject(embeddings=as.matrix(aligned_coor), key='aligned_')

### then, we perform UMAP or clustering based on the corrected PCA coordinates
obj = RunUMAP(object = obj, reduction = "aligned", dims = 1:30, min.dist = 0.3, verbose = FALSE)
obj = FindNeighbors(object = obj, reduction = "aligned", dims = 1:30, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 1, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 2, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 3, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 5, verbose = FALSE)
saveRDS(obj, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))




###############################################################
### Step-3: Annotating cell types based on the align_cds output

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))

anno = rep(NA, nrow(pd))
anno[pd$RNA_snn_res.1 %in% c(12)] = "Primordial germ cells"
anno[pd$RNA_snn_res.1 %in% c(6,16,21)] = "Epiblast"
anno[pd$RNA_snn_res.1 %in% c(20,23)] = "Neuroectoderm"
anno[pd$RNA_snn_res.1 %in% c(10)] = "Surface ectoderm"
anno[pd$RNA_snn_res.1 %in% c(7,15)] = "Endoderm"
anno[pd$RNA_snn_res.1 %in% c(0,2,4,5,11)] = "Lateral plate mesoderm"
anno[pd$RNA_snn_res.1 %in% c(3,8)] = "Cardiomyocytes"
anno[pd$RNA_snn_res.1 %in% c(1)] = "Paraxial mesoderm" ### Meox1, Meox2, Tbx1, Pax3
anno[pd$RNA_snn_res.1 %in% c(9,17)] = "Endothelial cells"
anno[pd$RNA_snn_res.1 %in% c(13,14)] = "Erythroid cells"
anno[pd$RNA_snn_res.1 %in% c(18,22)] = "Blood progenitors"
anno[pd$RNA_snn_res.1 %in% c(19)] = "Primitive streak"

pd$celltype = as.vector(anno)
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))



try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.12, color = "black") +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.08) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/plot/monoclonal_EB_TF_screen_celltype_aligned.png"),
               dpi = 300,
               height  = 6, 
               width = 6), silent = TRUE)


try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.12, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.08) +
        ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/plot/monoclonal_EB_TF_screen_celltype_label_aligned.png"),
               dpi = 300,
               height  = 6, 
               width = 6), silent = TRUE)

