
###################################################################################
### Integration with other datasets (mEB, natural embryo, or other in vitro models)
### Chengxiang Qiu
### Feb-20, 2025

###############################################
### Step-1: Integrating with wildtype mouse EBs

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

count = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/backup/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/obj_processed_pd.rds"))
count = count[,rownames(pd)]

### integrating with mEB_time_course_GEx
count_2 = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/count.rds"))
pd_2 = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/obj_sctransform_pd.rds"))
pd_2 = pd_2[,c("cell_id","experiment_id","celltype")]
pd_2$group = "mEB_time_course_GEx"

count_combine = cbind(count[rownames(count_2),], count_2)
pd_1 = pd[,c("cell_id","experiment_id")]
pd_1$celltype = pd$celltype
pd_1$group = "monoclonal_EB_proof"
pd_combine = rbind(pd_1, pd_2)

obj = CreateSeuratObject(count_combine, meta.data = pd_combine)
obj_processed = doClusterSeurat(obj)
saveRDS(obj_processed, paste0(work_path, "/analysis/monoclonal_EB_proof/obj_integration_mEB_time_course_GEx.rds"))


### making plots

pd = data.frame(obj_processed[[]])
pd$UMAP_1 = Embeddings(obj_processed, reduction = "umap")[,1]
pd$UMAP_2 = Embeddings(obj_processed, reduction = "umap")[,2]
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_proof/obj_integration_pd_mEB_time_course_GEx.rds"))

pd_x = pd %>% filter(group == "mEB_time_course_GEx") %>%
    group_by(celltype) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))
EB_celltype_color_code_sub = EB_celltype_color_code[names(EB_celltype_color_code) %in% pd_x$celltype]
EB_celltype_color_code_sub = data.frame(celltype = names(EB_celltype_color_code_sub), celltype_id = 1:length(EB_celltype_color_code_sub))
pd_x = pd_x %>% left_join(EB_celltype_color_code_sub, by = "celltype")

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "mEB_time_course_GEx"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "mEB_time_course_GEx"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype_id), color = "black", size = 6, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/obj_integration_mEB_time_course_GEx_1.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

pd_x = pd %>% filter(group == "monoclonal_EB_proof") %>%
    group_by(celltype) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))
EB_celltype_color_code_sub = EB_celltype_color_code[names(EB_celltype_color_code) %in% pd_x$celltype]
EB_celltype_color_code_sub = data.frame(celltype = names(EB_celltype_color_code_sub), celltype_id = 1:length(EB_celltype_color_code_sub))
pd_x = pd_x %>% left_join(EB_celltype_color_code_sub, by = "celltype")

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "monoclonal_EB_proof"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "monoclonal_EB_proof"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype_id), color = "black", size = 6, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/obj_integration_mEB_time_course_GEx_2.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)









############################################################
### Step-2: integration with Pijuan-Sala's data and JAX data


library(future)
library(future.apply)
plan("multicore", workers = 4)
options(future.globals.maxSize = 80000 * 1024^2)

### 1st dataset - Sam's EB data (time course)
count_1 = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/backup/count.rds"))
pd_1 = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/obj_processed_pd.rds"))
count_1 = count_1[,rownames(pd_1)]
pd_1$group = "monoclonal_EB_proof"


### 2nd dataset - our JAX data
work_path_2 = "/net/shendure/vol2/projects/cxqiu/work/jax/rna_seq"
pd_all = readRDS(paste0(work_path_2, "/mtx/adata_scale.obs.rds"))
pd = pd_all[pd_all$group %in% c("E85","E9","E95"),]
rownames(pd) = as.vector(pd$cell_id)
pd = pd[sample(1:nrow(pd), 200000),]
embryo_list = as.vector(names(table(pd$embryo_id)))

count_2 = NULL
for(j in embryo_list){
    
    print(paste0(j,"/",length(embryo_list)))
    count_j = readRDS(paste0(work_path_2, "/embryo/", j, "_gene_count.rds"))
    keep = colnames(count_j) %in% rownames(pd)
    count_2 = cbind(count_2, count_j[,keep])
    rm(count_j)
}
pd_2 = pd[colnames(count_2),]
pd_2$group = "jax"

fd = mouse_gene[rownames(count_2),]
fd_sub = fd[(fd$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & fd$chr %in% paste0("chr", c(1:19, "M")),]
count_2 = count_2[rownames(count_2) %in% as.vector(fd_sub$gene_ID),]



### 3rd dataset - Pijuan-Sala data
obj_3 = readRDS("/net/shendure/vol8/projects/cxqiu/data/mouse_gastrulation_Pijuan/Pijuan.rds")
pd_3 = data.frame(obj_3[[]])
pd_3$group = "pijuan"
count_3 = GetAssayData(obj_3, slot = "counts")


### merge three datasets
gene_overlap = intersect(rownames(count_1), intersect(rownames(count_2), rownames(count_3)))
print(length(gene_overlap))

obj_1 = CreateSeuratObject(count_1[gene_overlap,], meta.data = pd_1)
obj_2 = CreateSeuratObject(count_2[gene_overlap,], meta.data = pd_2)
obj_3 = CreateSeuratObject(count_3[gene_overlap,], meta.data = pd_3)

obj = merge(x = obj_1, y = c(obj_2, obj_3))

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

saveRDS(obj.integrated, paste0(work_path, "/analysis/monoclonal_EB_proof/obj_integration_embryo.rds"))

obj.integrated$UMAP_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
obj.integrated$UMAP_2 = Embeddings(obj.integrated, reduction = "umap")[,2]
saveRDS(data.frame(obj.integrated[[]]), paste0(work_path, "/analysis/monoclonal_EB_proof/obj_integration_pd_embryo.rds"))


### making plots 

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/obj_integration_pd_embryo.rds"))
pd_x = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/obj_processed_pd.rds"))

pd_1 = pd[pd$group == "jax",]
pd_1$celltypeannotation = as.vector(pd_1$major_trajectory)
pd_2 = pd[pd$group == "pijuan",]
pd_2$celltypeannotation = as.vector(pd_2$pre_celltype)
pd_3 = pd[pd$group == "monoclonal_EB_proof",]
pd_x = pd_x[rownames(pd_3),]
pd_3$celltypeannotation = as.vector(pd_x$celltype)

pd = rbind(pd_1, pd_2, pd_3)
pd = pd[,c("group","celltypeannotation", "UMAP_1", "UMAP_2")]
saveRDS(pd, paste0(work_path, "/analysis/monoclonal_EB_proof/obj_integration_pd_embryo.rds"))




pd_x = pd %>% filter(group == "monoclonal_EB_proof") %>%
    group_by(celltypeannotation) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))
EB_celltype_color_code_sub = EB_celltype_color_code[names(EB_celltype_color_code) %in% pd_x$celltypeannotation]
EB_celltype_color_code_sub = data.frame(celltypeannotation = names(EB_celltype_color_code_sub), celltype_id = 1:length(EB_celltype_color_code_sub))
pd_x = pd_x %>% left_join(EB_celltype_color_code_sub, by = "celltypeannotation")

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "monoclonal_EB_proof"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "monoclonal_EB_proof"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltypeannotation), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype_id), color = "black", size = 6, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/Integration_embryo_1.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

pd_x = pd %>% filter(group == "pijuan") %>%
    group_by(celltypeannotation) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))
pijuan_celltype_color_code_sub = pijuan_celltype_color_code[names(pijuan_celltype_color_code) %in% pd_x$celltypeannotation]
pijuan_celltype_color_code_sub = data.frame(celltypeannotation = names(pijuan_celltype_color_code_sub), celltype_id = 1:length(pijuan_celltype_color_code_sub))
pd_x = pd_x %>% left_join(pijuan_celltype_color_code_sub, by = "celltypeannotation")

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "pijuan"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "pijuan"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltypeannotation), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype_id), color = "black", size = 6, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=pijuan_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/Integration_embryo_2.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


pd_x = pd %>% filter(group == "jax") %>%
    group_by(celltypeannotation) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))
jax_celltype_color_code_sub = jax_celltype_color_code[names(jax_celltype_color_code) %in% pd_x$celltypeannotation]
jax_celltype_color_code_sub = data.frame(celltypeannotation = names(jax_celltype_color_code_sub), celltype_id = 1:length(jax_celltype_color_code_sub))
pd_x = pd_x %>% left_join(jax_celltype_color_code_sub, by = "celltypeannotation")


try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "jax"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "jax"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltypeannotation), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype_id), color = "black", size = 6, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=jax_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/Integration_embryo_3.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)









########################################################
### Step-3: Integration with Liberali gastruloid dataset

library(future)
library(future.apply)
plan("multicore", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2)

### 1st dataset - Sam's EB data (time course)
count_1 = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/backup/count.rds"))
pd_1 = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/obj_processed_pd.rds"))
count_1 = count_1[,rownames(pd_1)]
pd_1$group = "monoclonal_EB_proof"


mouse_gene <- read.table("~/work/tome/code/mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count_1 = count_1[rownames(count_1) %in% as.vector(mouse_gene_sub$gene_ID),]
fd_1 = mouse_gene[rownames(count_1),]
fd_1$rowSums = Matrix::rowSums(count_1)
fd_1_x = fd_1 %>% group_by(gene_short_name) %>% slice_max(order_by = rowSums, n = 1, with_ties = F)
count_1 = count_1[as.vector(fd_1_x$gene_ID),]
rownames(count_1) = as.vector(fd_1_x$gene_short_name)

### 2nd dataset - Liberali's data
count_2 = Matrix::readMM(paste0(data_path, "/GSE229513_UMI_counts.mtx"))
rownames(count_2) = as.vector(read.table(paste0(data_path, "/GSE229513_genes.tsv"), header=T)$gene)
colnames(count_2) = as.vector(read.table(paste0(data_path, "/GSE229513_barcodes.tsv"), sep="\t")$V1)
pd_2 = readRDS(paste0(data_path, "/df_cell.rds"))
pd_2 = pd_2[,c(1:4, 9, 13)]
pd_2$group = "Liberali"

### merge two datasets
gene_overlap = intersect(rownames(count_1), rownames(count_2))
print(length(gene_overlap))

obj_1 = CreateSeuratObject(count_1[gene_overlap,], meta.data = pd_1)
obj_2 = CreateSeuratObject(count_2[gene_overlap,], meta.data = pd_2)

obj = merge(x = obj_1, y = obj_2)

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

saveRDS(obj.integrated, paste0(work_path, "/analysis/monoclonal_EB_proof/obj_integration_Liberali.rds"))

obj.integrated$UMAP_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
obj.integrated$UMAP_2 = Embeddings(obj.integrated, reduction = "umap")[,2]
saveRDS(data.frame(obj.integrated[[]]), paste0(work_path, "/analysis/monoclonal_EB_proof/obj_integration_pd_Liberali.rds"))

### making plots

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/obj_integration_pd_Liberali.rds"))


try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "monoclonal_EB_proof"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "monoclonal_EB_proof"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype_id), color = "black", size = 6, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/obj_integration_mEB_time_course_GEx_2.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


pd_x = pd %>% filter(group == "Liberali") %>%
    group_by(celltypeannotation) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))
Liberali_celltype_color_code_sub = Liberali_celltype_color_code[names(Liberali_celltype_color_code) %in% pd_x$celltypeannotation]
Liberali_celltype_color_code_sub = data.frame(celltypeannotation = names(Liberali_celltype_color_code_sub), celltype_id = 1:length(Liberali_celltype_color_code_sub))
pd_x = pd_x %>% left_join(Liberali_celltype_color_code_sub, by = "celltypeannotation")

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "Liberali"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "Liberali"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltypeannotation), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype_id), color = "black", size = 6, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=Liberali_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/Integration_Liberali_1.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

pd_x = pd %>% filter(group == "monoclonal_EB_proof") %>%
    group_by(celltype) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))
EB_celltype_color_code_sub = EB_celltype_color_code[names(EB_celltype_color_code) %in% pd_x$celltype]
EB_celltype_color_code_sub = data.frame(celltype = names(EB_celltype_color_code_sub), celltype_id = 1:length(EB_celltype_color_code_sub))
pd_x = pd_x %>% left_join(EB_celltype_color_code_sub, by = "celltype")

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "monoclonal_EB_proof"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "monoclonal_EB_proof"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.1) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype_id), color = "black", size = 6, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/plot/Integration_Liberali_2.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)



