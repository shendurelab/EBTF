
######################################################
### Merging samples and performing dimension reduction
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]


###########################
### Step-1: merging samples

run_list = paste0("EB_", c(1:17,19,21:32))

obj = NULL
for(run_id in run_list){
    print(run_id)
    obj_tmp = readRDS(paste0(work_path, "/processing/all_TF_screen_GEx/", run_id, "/obj_all_TF_screen_GEx_", run_id, ".rds"))
    print(dim(obj_tmp))
    
    if(is.null(obj)){
        obj = obj_tmp
    } else {
        obj = merge(obj, obj_tmp)
    }
}

count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])

saveRDS(count, paste0(work_path, "/analysis/all_TF_screen_GEx/count.rds"))
saveRDS(pd, paste0(work_path, "/analysis/all_TF_screen_GEx/pd.rds"))


######################################################################
### Step-2: performing regular seurat analysis but including align_cds

count = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/backup/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/backup/pd.rds"))

### In the first try, I found there are two cell clusters showing lower UMI and higher MT%
pd_x = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/backup/obj_aligned_pd_first_try.rds"))
keep = !pd_x$RNA_snn_res.1 %in% c(17,26)

pd = pd[keep,]
count = count[,keep]

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
obj = FindClusters(object = obj, resolution = 2, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 3, verbose = FALSE)
saveRDS(obj, paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned_pd.rds"))




#################################
### Step-4: cell type annotations

anno = rep(NA, nrow(pd))
anno[pd$RNA_snn_res.1 %in% c(0,5,8,9,12)] = "Primordial germ cells"
anno[pd$RNA_snn_res.1 %in% c(4,10)] = "Epiblast"
anno[pd$RNA_snn_res.1 %in% c(14)] = "Primitive streak"
anno[pd$RNA_snn_res.1 %in% c(28)] = "ESCs (2-cell state)"
anno[pd$RNA_snn_res.1 %in% c(21)] = "Neuroectoderm"
anno[pd$RNA_snn_res.1 %in% c(18)] = "Surface ectoderm"
anno[pd$RNA_snn_res.1 %in% c(24)] = "Notochord"
anno[pd$RNA_snn_res.1 %in% c(15)] = "Definitive endoderm"
anno[pd$RNA_snn_res.1 %in% c(11,16)] = "Extraembryonic visceral endoderm"
anno[pd$RNA_snn_res.1 %in% c(1,2,3,6,7,17,22,26,29)] = "Parietal endoderm"
anno[pd$RNA_snn_res.1 %in% c(13,19,20,27)] = "Lateral plate mesoderm"
anno[pd$RNA_snn_res.1 %in% c(23)] = "Endothelial cells"
anno[pd$RNA_snn_res.1 %in% c(25)] = "Erythroid cells"

anno[pd$RNA_snn_res.2 %in% c(4,16,17)] = "Epiblast"
anno[pd$RNA_snn_res.2 %in% c(3,7,8,9,18,26)] = "Primordial germ cells"

anno[pd$RNA_snn_res.3 %in% c(56)] = "Cardiomyocytes"
anno[pd$RNA_snn_res.3 %in% c(49)] = "Nascent mesoderm"
anno[pd$RNA_snn_res.3 %in% c(55)] = "Blood progenitors"

pd$celltype = as.vector(anno)
saveRDS(pd, paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned_pd.rds"))


try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.12, color = "black") +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.08) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/all_TF_screen_GEx_celltype.png"),
               dpi = 300,
               height  = 7, 
               width = 8), silent = TRUE)


try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.12, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.08) +
        ggrepel::geom_text_repel(data = pd %>% group_by(celltype) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = celltype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/all_TF_screen_GEx_celltype_label.png"),
               dpi = 300,
               height  = 7, 
               width = 8), silent = TRUE)




