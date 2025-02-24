
#####################################################
### Merging three timepoints and performing embedding
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene <- read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

####################################################
### Step-1: integration with natural embryo datasets

library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 80000 * 1024^2)

### 1st dataset - our EB data (time course)
count_1 = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/count.rds"))
pd_1 = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/pd.rds"))
pd_1$group = "mEB_time_course_GEx"

### 2nd dataset - natural embryos during organogenesis
### reference: https://pubmed.ncbi.nlm.nih.gov/38355799/
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

### 3rd dataset - natural embryos during gastrulation
### reference: https://pubmed.ncbi.nlm.nih.gov/30787436/
obj_3 = readRDS("./Pijuan.rds")
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

saveRDS(obj.integrated, paste0(work_path, "/analysis/mEB_time_course_GEx/obj_integration_embryo.rds"))

obj.integrated$UMAP_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
obj.integrated$UMAP_2 = Embeddings(obj.integrated, reduction = "umap")[,2]
saveRDS(data.frame(obj.integrated[[]]), paste0(work_path, "/analysis/mEB_time_course_GEx/obj_integration_pd_embryo.rds"))


##############################################################
### Step-2: making three UMAP plots for the integration result

pd = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/obj_integration_pd_embryo.rds"))
pd_x = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/obj_sctransform_pd.rds"))

pd_1 = pd[pd$group == "jax",]
pd_1$celltypeannotation = as.vector(pd_1$major_trajectory)
pd_2 = pd[pd$group == "pijuan",]
pd_2$celltypeannotation = as.vector(pd_2$pre_celltype)
pd_3 = pd[pd$group == "mEB_time_course_GEx",]
pd_x = pd_x[rownames(pd_3),]
pd_3$celltypeannotation = as.vector(pd_x$celltype)

pd = rbind(pd_1, pd_2, pd_3)
pd = pd[,c("group","celltypeannotation", "UMAP_1", "UMAP_2")]
saveRDS(pd, paste0(work_path, "/analysis/mEB_time_course_GEx/obj_integration_pd_embryo.rds"))

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "mEB_time_course_GEx"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "mEB_time_course_GEx"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltypeannotation), size=0.1) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/Integration_embryo_1.png"),
               dpi = 300,
               height  = 6, 
               width = 6), silent = TRUE)


try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "pijuan"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "pijuan"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltypeannotation), size=0.1) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=pijuan_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/Integration_embryo_2.png"),
               dpi = 300,
               height  = 6, 
               width = 6), silent = TRUE)


try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "jax"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "jax"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltypeannotation), size=0.1) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=jax_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/Integration_embryo_3.png"),
               dpi = 300,
               height  = 6, 
               width = 6), silent = TRUE)




######################################################################################
### Step-3: Performing NNLS to align cell types between mEBs and natural embryos (JAX)

pd_integration = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/obj_integration_pd_embryo.rds"))

pd_all = readRDS(paste0(work_path_2, "/mtx/adata_scale.obs.rds"))
pd = pd_all[rownames(pd_all) %in% rownames(pd_integration),]
embryo_list = as.vector(names(table(pd$embryo_id)))

count_2 = NULL
for(j in embryo_list){
    print(paste0(j,"/",length(embryo_list)))
    count_j = readRDS(paste0(work_path_2, "/embryo/", j, "_gene_count.rds"))
    keep = colnames(count_j) %in% rownames(pd)
    count_2 = cbind(count_2, count_j[,keep,drop=FALSE])
    rm(count_j)
}
pd_2 = pd[colnames(count_2),]

count_2 = t(t(count_2) / colSums(count_2)) * 10000
count_2@x = log(count_2@x + 1)

major_trajectory_list = unique(pd_2$major_trajectory)

exp_aggr = NULL
for(i in major_trajectory_list){
    exp_aggr = cbind(exp_aggr, Matrix::rowMeans(count_2[,pd_2$major_trajectory == i,drop=FALSE]))
}
colnames(exp_aggr) = major_trajectory_list

saveRDS(list(exp_aggr, pd_2), paste0(work_path, "/analysis/mEB_time_course_GEx/JAX_exp.rds"))

### performing NNLs followed by making heatmap

dat_jax = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/JAX_exp.rds"))
dat_jax = dat_jax[[1]]
dat_mEB = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/mEB_exp.rds"))

gene_overlap = intersect(rownames(dat_jax), rownames(dat_mEB))
dat_jax = dat_jax[gene_overlap,]
dat_mEB = dat_mEB[gene_overlap,]

conn = correlation_analysis_bidirection(as.matrix(dat_jax), as.matrix(dat_mEB), fold.change = 1.5, top_gene_num = 1500, spec_gene_num = 1500)
conn$beta = 2*(conn$beta_1+0.01)*(conn$beta_2+0.01)
saveRDS(conn, paste0(work_path, "/analysis/mEB_time_course_GEx/NNLS_jax_mEB.rds"))

dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]

# manually order rownames (cell types)
col_names = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% colnames(dat)]
row_names = NULL
for(i in col_names){
    for(j in 1:nrow(dat)){
        tmp = names(sort(dat[j,], decreasing=T))[1]
        if(tmp == i & !rownames(dat)[j] %in% row_names){
            row_names = c(row_names, rownames(dat)[j])
        }
    }
}
dat = dat[row_names,col_names]

pdf(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/NNLS_jax_mEB.pdf"),12,12)
heatmap.2(as.matrix(t(dat)), 
          col=viridis, 
          scale="row", 
          Rowv = FALSE, 
          Colv = FALSE, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(15,15))
dev.off()


##################
### making heatmap 

dat = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/JAX_exp.rds"))
exp_aggr = dat[[1]]
pd = dat[[2]]

res = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/celltype_TFs.rds"))

res_top_3_fc = res %>% group_by(cluster) %>% slice_max(order_by = avg_logFC, n = 3, with_ties = FALSE) %>% 
    select(cluster, gene_short_name, gene_ID) %>% as.data.frame()

EB_celltype_list = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% res_top_3_fc$cluster]

res_top_3 = res_top_3_fc
res_top_3$cluster = factor(res_top_3$cluster, levels = EB_celltype_list)
res_top_3 = res_top_3[order(res_top_3$cluster),]

res_top_3_gene = unique(res_top_3[,c("gene_short_name", "gene_ID")])
exp_aggr = exp_aggr[as.vector(res_top_3_gene$gene_ID),]
rownames(exp_aggr) = as.vector(res_top_3_gene$gene_short_name)

col_names = c("Neuroectoderm_and_glia", "CNS_neurons", "Epithelium", "Ependymal_cells", 
              "Intestine", "Mesoderm", "Cardiomyocytes", "Endothelium", "White_blood_cells", "Primitive_erythroid")

exp_aggr = exp_aggr[,col_names]

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/celltype_TFs_JAX_major_trajectory.pdf"), 8, 5)
heatmap.2(as.matrix(t(exp_aggr)), 
          col=Colors, 
          scale="col", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()

write.table(rownames(exp_aggr), paste0(work_path, "/analysis/mEB_time_course_GEx/plot/celltype_TFs_JAX_major_trajectory.colnames.txt"), row.names=F, col.names=F, sep="\t", quote=F)





###############################################################################################
### Step-4: Performing NNLS to align cell types between mEBs and natural embryos (gastrulation)

obj_3 = readRDS("/net/shendure/vol10/projects/cxqiu/nobackup/data/mouse_gastrulation_Pijuan/Pijuan.rds")
pd_3 = data.frame(obj_3[[]])
pd_3$group = "pijuan"
count_3 = GetAssayData(obj_3, slot = "counts")

count_3 = t(t(count_3) / colSums(count_3)) * 10000
count_3@x = log(count_3@x + 1)

pre_celltype_list = unique(pd_3$pre_celltype)

exp_aggr = NULL
for(i in pre_celltype_list){
    exp_aggr = cbind(exp_aggr, Matrix::rowMeans(count_3[,pd_3$pre_celltype == i,drop=FALSE]))
}
colnames(exp_aggr) = pre_celltype_list

saveRDS(list(exp_aggr, pd_3), paste0(work_path, "/analysis/mEB_time_course_GEx/Pijuan_exp.rds"))


##############################################
### performing NNLs followed by making heatmap

dat_pijuan = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/Pijuan_exp.rds"))
dat_pijuan = dat_pijuan[[1]]
dat_mEB = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/mEB_exp.rds"))

gene_overlap = intersect(rownames(dat_pijuan), rownames(dat_mEB))
dat_pijuan = dat_pijuan[gene_overlap,]
dat_mEB = dat_mEB[gene_overlap,]

conn = correlation_analysis_bidirection(as.matrix(dat_pijuan), as.matrix(dat_mEB), fold.change = 1.5, top_gene_num = 1500, spec_gene_num = 1500)
conn$beta = 2*(conn$beta_1+0.01)*(conn$beta_2+0.01)
saveRDS(conn, paste0(work_path, "/analysis/mEB_time_course_GEx/NNLS_pijuan_mEB.rds"))

dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]

# manually order rownames (cell types)
col_names = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% colnames(dat)]
row_names = NULL
for(i in col_names){
    for(j in 1:nrow(dat)){
        tmp = names(sort(dat[j,], decreasing=T))[1]
        if(tmp == i & !rownames(dat)[j] %in% row_names){
            row_names = c(row_names, rownames(dat)[j])
        }
    }
}
dat = dat[row_names,col_names]

pdf(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/NNLS_pijuan_mEB.pdf"),12,12)
heatmap.2(as.matrix(t(dat)), 
          col=viridis, 
          scale="row", 
          Rowv = FALSE, 
          Colv = FALSE, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(15,15))
dev.off()





########################################################
### Step-5: Integration with Liberali gastruloid dataset


library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 30000 * 1024^2)

### 1st dataset - Our EB data (time course)
count_1 = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/count.rds"))
pd_1 = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/pd.rds"))
pd_1$group = "mEB_time_course_GEx"

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count_1 = count_1[rownames(count_1) %in% as.vector(mouse_gene_sub$gene_ID),]
fd_1 = mouse_gene[rownames(count_1),]
fd_1$rowSums = Matrix::rowSums(count_1)
fd_1_x = fd_1 %>% group_by(gene_short_name) %>% slice_max(order_by = rowSums, n = 1, with_ties = F)
count_1 = count_1[as.vector(fd_1_x$gene_ID),]
rownames(count_1) = as.vector(fd_1_x$gene_short_name)

### 2nd dataset - Suppinger's dataset
### reference: https://pubmed.ncbi.nlm.nih.gov/37209681/
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

saveRDS(obj.integrated, paste0(work_path, "/analysis/mEB_time_course_GEx/obj_integration_Liberali.rds"))

obj.integrated$UMAP_1 = Embeddings(obj.integrated, reduction = "umap")[,1]
obj.integrated$UMAP_2 = Embeddings(obj.integrated, reduction = "umap")[,2]
saveRDS(data.frame(obj.integrated[[]]), paste0(work_path, "/analysis/mEB_time_course_GEx/obj_integration_pd_Liberali.rds"))


################
### making plots

pd = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/obj_integration_pd_Liberali.rds"))

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "Liberali"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "Liberali"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltypeannotation), size=0.1) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=Liberali_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/Integration_Liberali_1.png"),
               dpi = 300,
               height  = 6, 
               width = 6), silent = TRUE)

try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(pd, group == "mEB_time_course_GEx"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
        geom_point(data = subset(pd, group == "mEB_time_course_GEx"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltypeannotation), size=0.1) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/Integration_Liberali_2.png"),
               dpi = 300,
               height  = 6, 
               width = 6), silent = TRUE)



