
######################################################
### Merging samples and performing dimension reduction
### Chengxiang Qiu
### Feb-20, 2025

###########################
### Step-1: merging samples

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

run_list = paste0("EB_", c(1:32))

obj = NULL
for(run_id in run_list){
    print(run_id)
    obj_tmp = readRDS(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/obj_EB21_CRISPRcut_125TF_GEx_", run_id, ".rds"))
    print(dim(obj_tmp))
    
    if(is.null(obj)){
        obj = obj_tmp
    } else {
        obj = merge(obj, obj_tmp)
    }
}

count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])

saveRDS(count, paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/backup/count.rds"))
saveRDS(pd, paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/backup/pd.rds"))


######################################################################
### Step-2: performing regular seurat analysis but including align_cds

count = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/backup/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/backup/pd.rds"))

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
saveRDS(obj, paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/backup/obj_aligned_first_try.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/backup/obj_aligned_pd_first_try.rds"))




#####################################################
### Step-3: excluding some clusters with very low UMI

count = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/backup/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/backup/pd.rds"))

### In the first try, I found there are several cell clusters showing lower UMI 
pd_x = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/backup/obj_aligned_pd_first_try.rds"))
keep = !pd_x$RNA_snn_res.1 %in% c(0,9,24,25,27,28,29,30,31,32)

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
obj = FindClusters(object = obj, resolution = 5, verbose = FALSE)
saveRDS(obj, paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned_pd.rds"))



#################################
### Step-4: cell type annotations

anno = rep(NA, nrow(pd))
anno[pd$RNA_snn_res.1 %in% c(0,1,2,3,6,7,12,17,18)] = "Primordial germ cells"
anno[pd$RNA_snn_res.1 %in% c(4,5)] = "Epiblast"
anno[pd$RNA_snn_res.1 %in% c(19)] = "Primitive streak"
anno[pd$RNA_snn_res.1 %in% c(22)] = "ESCs (2-cell state)"
anno[pd$RNA_snn_res.1 %in% c(8)] = "Neuroectoderm"
anno[pd$RNA_snn_res.1 %in% c(23)] = "Eye field"
anno[pd$RNA_snn_res.1 %in% c(15)] = "Early neurons"
anno[pd$RNA_snn_res.1 %in% c(20)] = "Floor plate"
anno[pd$RNA_snn_res.1 %in% c(11)] = "Surface ectoderm"
anno[pd$RNA_snn_res.1 %in% c(16)] = "Definitive endoderm"
anno[pd$RNA_snn_res.1 %in% c(14)] = "Extraembryonic visceral endoderm"
anno[pd$RNA_snn_res.1 %in% c(9,13,27)] = "Parietal endoderm"
anno[pd$RNA_snn_res.1 %in% c(10,21,24)] = "Lateral plate mesoderm"
anno[pd$RNA_snn_res.1 %in% c(25)] = "Endothelial cells"
anno[pd$RNA_snn_res.1 %in% c(26)] = "Erythroid cells"

anno[pd$RNA_snn_res.1 %in% c(26) & pd$RNA_snn_res.5 %in% c(80)] = "Blood progenitors"
anno[pd$RNA_snn_res.1 %in% c(16) & pd$RNA_snn_res.3 %in% c(50)] = "Extraembryonic visceral endoderm"
anno[pd$RNA_snn_res.1 %in% c(21) & pd$RNA_snn_res.3 %in% c(58)] = "Cardiomyocytes"
anno[pd$RNA_snn_res.1 %in% c(21) & pd$RNA_snn_res.3 %in% c(52)] = "Nascent mesoderm"
anno[pd$RNA_snn_res.1 %in% c(19) & pd$RNA_snn_res.5 %in% c(69)] = "Notochord"

pd$celltype = as.vector(anno)
saveRDS(pd, paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned_pd.rds"))


try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.12, color = "black") +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.08) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/EB21_CRISPRcut_125TF_GEx_celltype.png"),
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
        ggsave(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/EB21_CRISPRcut_125TF_GEx_celltype_label.png"),
               dpi = 300,
               height  = 7, 
               width = 8), silent = TRUE)

pd$replicate = rep("rep1", nrow(pd))
pd$replicate[pd$experiment_id %in% paste0("EB21_CRISPRcut_125TF_GEx_EB_", c(17:32))] = "rep2"

rep_color = c("rep1"="red", "rep2"="blue")
for(i in c("rep1", "rep2")){
    try(ggplot() +
            geom_point(data = subset(pd, replicate != i), aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
            geom_point(data = subset(pd, replicate == i), aes(x = UMAP_1, y = UMAP_2), size=0.2, color = rep_color[i]) +
            theme_void() +
            theme(legend.position="none") +
            ggsave(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/EB21_CRISPRcut_125TF_GEx_celltype_", i, ".png"),
                   dpi = 300,
                   height  = 7, 
                   width = 8), silent = TRUE)
}









