
###########################################################################################
### Performing DEGs between Hand1 and NTCs at pseudobulk level, in the "arrayed" experiment
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]


#################################################################
### Step-1: Assigning gRNAs and targets for individual colonotype

cell_gRNA = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep1.rds"))
cell_gRNA$target = unlist(lapply(as.vector(cell_gRNA$gRNA), function(x) strsplit(x,"[_]")[[1]][1])) 

cell_gRNA_uniq = cell_gRNA %>% filter(pct >= 70) %>% arrange(clonotype_id)
cell_gRNA_mult = cell_gRNA %>% filter(!clonotype_id %in% as.vector(cell_gRNA_uniq$clonotype_id), order %in% c(1,2))

### out of 202 clonotypes, 187 are assigned with uique gRNAs, 15 are assigned with double gRNAs
cell_gRNA = rbind(cell_gRNA_uniq, cell_gRNA_mult)
cell_target = cell_gRNA %>% select(clonotype_id, target) %>% unique()

### we further set a new cutoff of cell number for identifying EBs (cell # = 100)
cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep1.rds"))
cell_BC_assigned_num = cell_BC_assigned %>% group_by(clonotype_id) %>% tally() %>% filter(n >= 100)
cell_target = cell_target[cell_target$clonotype_id %in% as.vector(cell_BC_assigned_num$clonotype_id),]
cell_BC_assigned = cell_BC_assigned[cell_BC_assigned$clonotype_id %in% as.vector(cell_BC_assigned_num$clonotype_id),]


########################################
### Step-2: creating pseudobulk profiles

obj = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))
count = GetAssayData(obj, slot = "counts")

clonotype_list = unique(cell_target$clonotype_id)
count_aggr = NULL
for(i in 1:length(clonotype_list)){
    print(paste0(i,"/",length(clonotype_list)))
    clonotype_i = clonotype_list[i]
    count_sub = count[,colnames(count) %in% as.vector(cell_BC_assigned$cell_id[cell_BC_assigned$clonotype_id == clonotype_i])]
    count_aggr = cbind(count_aggr, Matrix::rowSums(count_sub))
}
count_aggr = as.matrix(count_aggr)
colnames(count_aggr) = clonotype_list
saveRDS(list(count_aggr, cell_target), paste0(work_path, "/analysis/monoclonal_EB_TF_screen/DEG/count_aggr_1.rds"))


######################################
### Step-3: Performing DESeq2 analysis

library(DESeq2)
source("~/work/scripts/utils.R")
mouse_gene <- read.table("~/work/tome/code/mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf"

dat_tmp = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/DEG/count_aggr_1.rds"))
count_aggr = dat_tmp[[1]]
cell_target = dat_tmp[[2]]
rownames(cell_target) = as.vector(cell_target$clonotype_id)
print(sum(colnames(count_aggr) != rownames(cell_target)))

########
target = "Hand1"

keep = cell_target$target %in% c(target, "NONTARGETING")
metadata = cell_target[keep,]
count_data = count_aggr[,keep]
metadata$condition = if_else(metadata$target == target, "KO", "Control")
metadata$condition <- factor(metadata$condition, levels = c("Control", "KO"))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata,
                              design = ~ condition)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
results <- results(dds, contrast = c("condition", "KO", "Control"))
summary(results)

results = as.data.frame(results)
results_sig =  results %>% mutate(gene_ID = rownames(results)) %>% 
    filter(!is.na(padj), padj < 0.05) %>%
    left_join(mouse_gene[,c("gene_ID","gene_short_name")], by = "gene_ID")

write.table(results_sig, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/DEG/Hand1_vs_NTC_1.txt"), row.names=F, quote=F, sep="\t")



##########################################################################
### Step-4: performing a single-cell level analysis (w/o monoclonal types)

obj_orig = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))

gRNA_matrix = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/gRNA_matrix_rep1.rds"))
gRNA_matrix = gRNA_matrix[,colSums(gRNA_matrix) >= 10]
dat = t(t(gRNA_matrix)/colSums(gRNA_matrix))
target_list = unlist(lapply(rownames(dat), function(x) strsplit(x,"[_]")[[1]][1]))
Hand1_cell_list = colnames(dat)[colSums(dat[target_list == "Hand1",]) >= 0.9]
NTC_cell_list = colnames(dat)[colSums(dat[target_list == "NONTARGETING",]) >= 0.9]

count = GetAssayData(obj_orig, slot = "counts")
count = count[,colnames(count) %in% c(Hand1_cell_list, NTC_cell_list)]
count = count[rowSums(count) > 10,]
obj = CreateSeuratObject(count)
Idents(obj) = if_else(colnames(obj) %in% Hand1_cell_list, "KO", "Control")
obj = NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

res = FindMarkers(obj, ident.1 = "KO", ident.2 = "Control")
res_sig = res[res$p_val_adj < 0.05,]
x_sc_up = rownames(res_sig)[res_sig$avg_log2FC > 0]
x_sc_down = rownames(res_sig)[res_sig$avg_log2FC < 0]

result_sig = read.table(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/DEG/Hand1_vs_NTC_1.txt"), header=T)
x_bulk_up = as.vector(result_sig$gene_ID[result_sig$log2FoldChange > 0])
x_bulk_down = as.vector(result_sig$gene_ID[result_sig$log2FoldChange < 0])

library(ggvenn)
library(patchwork)
p1 = ggvenn(list(sc_up = x_sc_up, bulk_up = x_bulk_up), show_percentage = T)
p2 = ggvenn(list(sc_down = x_sc_down, bulk_down = x_bulk_down), show_percentage = T)
ggsave("~/share/EBTF_experiment1_Venn.pdf", p1+p2, width = 8, height = 4)


################################################################################
### Step-5: performing GO analysis on shared genes vs. only found in single-cell

library(topGO)

res = NULL

overlap_up = x_sc_up[x_sc_up %in% x_bulk_up]
overlap_up_name = as.vector(mouse_gene_sub$gene_short_name[mouse_gene_sub$gene_ID %in% overlap_up])
geneList = ifelse(unique(as.vector(mouse_gene_sub$gene_short_name)) %in% overlap_up_name, 1, 0)
names(geneList) = unique(as.vector(mouse_gene_sub$gene_short_name))
print(table(geneList))

GOdata <- new("topGOdata",
              ontology = "BP", 
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
res_i = GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 200)
res_i$gene_group = "Overlap_Up"
res = rbind(res, res_i)

sc_only_up = x_sc_up[!x_sc_up %in% x_bulk_up]
sc_only_up_name = as.vector(mouse_gene_sub$gene_short_name[mouse_gene_sub$gene_ID %in% sc_only_up])
geneList = ifelse(unique(as.vector(mouse_gene_sub$gene_short_name)) %in% sc_only_up_name, 1, 0)
names(geneList) = unique(as.vector(mouse_gene_sub$gene_short_name))
print(table(geneList))

GOdata <- new("topGOdata",
              ontology = "BP", 
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
res_i = GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 200)
res_i$gene_group = "SC_Only_Up"
res = rbind(res, res_i)

overlap_down = x_sc_down[x_sc_down %in% x_bulk_down]
overlap_down_name = as.vector(mouse_gene_sub$gene_short_name[mouse_gene_sub$gene_ID %in% overlap_down])
geneList = ifelse(unique(as.vector(mouse_gene_sub$gene_short_name)) %in% overlap_down_name, 1, 0)
names(geneList) = unique(as.vector(mouse_gene_sub$gene_short_name))
print(table(geneList))

GOdata <- new("topGOdata",
              ontology = "BP", 
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
res_i = GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 200)
res_i$gene_group = "Overlap_Down"
res = rbind(res, res_i)

sc_only_down = x_sc_down[!x_sc_down %in% x_bulk_down]
sc_only_down_name = as.vector(mouse_gene_sub$gene_short_name[mouse_gene_sub$gene_ID %in% sc_only_down])
geneList = ifelse(unique(as.vector(mouse_gene_sub$gene_short_name)) %in% sc_only_down_name, 1, 0)
names(geneList) = unique(as.vector(mouse_gene_sub$gene_short_name))
print(table(geneList))

GOdata <- new("topGOdata",
              ontology = "BP", 
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
res_i = GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 200)
res_i$gene_group = "SC_Only_Down"
res = rbind(res, res_i)

res$experiment = "experiment_1"

write.table(res, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/DEG/DEG_GO_analysis_1.txt"), row.names=F, sep="\t", quote=F)
