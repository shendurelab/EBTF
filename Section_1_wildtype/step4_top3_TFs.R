
###########################################################
### Comparing specific TFs between natural embryos and mEBs
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene <- read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

#############################################
### Step-1: Identifying top 3 TF in mouse EBs

run_id = "day7"
obj_1 = readRDS(paste0(paste0(work_path, "/processing/mEB_time_course_GEx/", run_id, "/obj_", "mEB_time_course_GEx_", run_id, ".rds")))
run_id = "day14"
obj_2 = readRDS(paste0(paste0(work_path, "/processing/mEB_time_course_GEx/", run_id, "/obj_", "mEB_time_course_GEx_", run_id, ".rds")))
run_id = "day21"
obj_3 = readRDS(paste0(paste0(work_path, "/processing/mEB_time_course_GEx/", run_id, "/obj_", "mEB_time_course_GEx_", run_id, ".rds")))

obj = merge(obj_1, c(obj_2, obj_3))
count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])

pd_sctransform = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/obj_sctransform_pd.rds"))
pd_sctransform = pd_sctransform[rownames(pd),]
pd$celltype = as.vector(pd_sctransform$celltype)

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]

obj = CreateSeuratObject(count, meta.data = pd)
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 10000)
HVF = VariableFeatures(obj)

mouse_gene_sub = subset(mouse_gene, gene_ID %in% rownames(obj))

tf_list = read.table("~/work/scripts/sam_tf/scaled_gRNA_target_name_protospacer_cut.txt", as.is=T)
names(tf_list) = c("gene_short_name", "barcode")
gene_include = as.vector(mouse_gene_sub$gene_ID[mouse_gene_sub$gene_short_name %in% as.vector(tf_list$gene_short_name)])

Idents(obj) = as.vector(obj$celltype)
obj_sub = subset(obj, downsample = 1000)

res = FindAllMarkers(obj_sub, features = gene_include, logfc.threshold = 0.1, min.pct = 0.1, only.pos = T)
res = res %>% rename(gene_ID = gene) %>% left_join(mouse_gene, by = "gene_ID") %>% filter(p_val_adj < 0.05) %>% as.data.frame()
saveRDS(res, paste0(work_path, "/analysis/mEB_time_course_GEx/celltype_TFs.rds"))

res_top_3_pval = res %>% group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 3, with_ties = FALSE) %>% 
    select(cluster, gene_short_name, gene_ID) %>% as.data.frame()

res_top_3_fc = res %>% group_by(cluster) %>% slice_max(order_by = avg_logFC, n = 3, with_ties = FALSE) %>% 
    select(cluster, gene_short_name, gene_ID) %>% as.data.frame()

EB_celltype_list = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% res_top_3$cluster]

res_top_3 = res_top_3_fc
res_top_3$cluster = factor(res_top_3$cluster, levels = EB_celltype_list)
res_top_3 = res_top_3[order(res_top_3$cluster),]

### making heatmap

exp = GetAssayData(obj, slot = "counts")
exp = t(t(exp) / colSums(exp)) * 10000
exp@x = log(exp@x + 1)

exp_aggr = NULL
for(i in EB_celltype_list){
    exp_aggr = cbind(exp_aggr, Matrix::rowMeans(exp[,obj$celltype == i]))
}
colnames(exp_aggr) = EB_celltype_list
saveRDS(exp_aggr, paste0(work_path, "/analysis/mEB_time_course_GEx/mEB_exp.rds"))

res_top_3_gene = unique(res_top_3[,c("gene_short_name", "gene_ID")])
exp_aggr = exp_aggr[as.vector(res_top_3_gene$gene_ID),]
rownames(exp_aggr) = as.vector(res_top_3_gene$gene_short_name)

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/celltype_TFs.pdf"), 8, 5)
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

write.table(rownames(exp_aggr), paste0(work_path, "/analysis/mEB_time_course_GEx/plot/celltype_TFs.colnames.txt"), row.names=F, col.names=F, sep="\t", quote=F)


####################################################
### Step-2: Plotting top 3 TF in Pijuan-Sala dataset

dat = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/Pijuan_exp.rds"))
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

### aligned cell types
conn = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/NNLS_pijuan_mEB.rds"))
conn_sub = conn %>% filter(beta > 0.1)
conn_sub$target = factor(conn_sub$target, levels = EB_celltype_list[EB_celltype_list %in% conn_sub$target])
conn_sub = conn_sub %>% arrange(target)

col_names = c("PGC","Epiblast", "Forebrain/Midbrain/Hindbrain", "Surface ectoderm", 
              "Notochord", "Def. endoderm", "ExE endoderm", "Parietal endoderm",
              "Mesenchyme", "Cardiomyocytes", "Endothelium", "Erythroid3")

exp_aggr = exp_aggr[,col_names]

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/celltype_TFs_Pijuan_pre_celltype.pdf"), 8, 5)
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

write.table(rownames(exp_aggr), paste0(work_path, "/analysis/mEB_time_course_GEx/plot/celltype_TFs_Pijuan_pre_celltype.colnames.txt"), row.names=F, col.names=F, sep="\t", quote=F)



################################################################################
### Step-3: Identify top 3 TFs in Pijuan-Sala data and then plot them in the mEB


obj_3 = readRDS("/net/shendure/vol8/projects/cxqiu/data/mouse_gastrulation_Pijuan/Pijuan.rds")
pd_3 = data.frame(obj_3[[]])
pd_3$group = "pijuan"
count_3 = GetAssayData(obj_3, slot = "counts")

col_names = c("PGC","Epiblast", "Forebrain/Midbrain/Hindbrain", "Surface ectoderm", 
              "Notochord", "Def. endoderm", "ExE endoderm", "Parietal endoderm",
              "Mesenchyme", "Cardiomyocytes", "Endothelium", "Erythroid3")

keep = pd_3$pre_celltype %in% col_names
obj = CreateSeuratObject(count_3[,keep], meta.data = pd_3[keep,])
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)


run_id = "day7"
obj_1 = readRDS(paste0(paste0(work_path, "/processing/mEB_time_course_GEx/", run_id, "/obj_", "mEB_time_course_GEx_", run_id, ".rds")))


mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
mouse_gene_sub = subset(mouse_gene_sub, gene_ID %in% rownames(obj_1))
tf_list = read.table("~/work/scripts/sam_tf/scaled_gRNA_target_name_protospacer_cut.txt", as.is=T)
names(tf_list) = c("gene_short_name", "barcode")
gene_include = as.vector(mouse_gene_sub$gene_ID[mouse_gene_sub$gene_short_name %in% as.vector(tf_list$gene_short_name)])
gene_include = gene_include[gene_include %in% rownames(obj)]

Idents(obj) = as.vector(obj$pre_celltype)
obj_sub = subset(obj, downsample = 1000)

res = FindAllMarkers(obj_sub, features = gene_include, logfc.threshold = 0.1, min.pct = 0.1, only.pos = T)
res = res %>% rename(gene_ID = gene) %>% left_join(mouse_gene, by = "gene_ID") %>% filter(p_val_adj < 0.05) %>% as.data.frame()

### plot the top3 TFs in each cell type in Pijuan-Sala'data

exp_aggr = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/Pijuan_exp.rds"))[[1]]

res_top_3_fc = res %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 3, with_ties = FALSE) %>% 
    select(cluster, gene_short_name, gene_ID) %>% as.data.frame()

res_top_3 = res_top_3_fc
res_top_3$cluster = factor(res_top_3$cluster, levels = col_names)
res_top_3 = res_top_3[order(res_top_3$cluster),]

res_top_3_gene = unique(res_top_3[,c("gene_short_name", "gene_ID")])
exp_aggr = exp_aggr[as.vector(res_top_3_gene$gene_ID),]
rownames(exp_aggr) = as.vector(res_top_3_gene$gene_short_name)

exp_aggr = exp_aggr[,col_names]

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/celltype_TFs_Pijuan_pre_celltype_opposite.pdf"), 8, 5)
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

write.table(rownames(exp_aggr), paste0(work_path, "/analysis/mEB_time_course_GEx/plot/celltype_TFs_Pijuan_pre_celltype_opposite.colnames.txt"), row.names=F, col.names=F, sep="\t", quote=F)

### plot the top3 TFs in each cell type in mEB dataset


exp_aggr = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/mEB_exp.rds"))

EB_celltype_list = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% colnames(exp_aggr)]

exp_aggr = exp_aggr[as.vector(res_top_3_gene$gene_ID),EB_celltype_list]
rownames(exp_aggr) = as.vector(res_top_3_gene$gene_short_name)

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/celltype_TFs_opposite.pdf"), 8, 5)
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





