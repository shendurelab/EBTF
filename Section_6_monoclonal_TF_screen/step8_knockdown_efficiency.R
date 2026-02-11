
##################################################################
### In the monoclone screen expression, checking the KO efficiency

### Here, we treated NTC and all the other targets as "control" group

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

pd_all = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
obj = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))
count_all = GetAssayData(obj, slot = "counts")


#############
### Arrayed

hooke_ccs = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/hooke_ccs_pd.rds"))[[1]]
hooke_ccs = data.frame(clonotype_id = rownames(hooke_ccs), target = hooke_ccs$target)
cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep1.rds"))
cell_BC_assigned = cell_BC_assigned[cell_BC_assigned$clonotype_id %in% hooke_ccs$clonotype_id,]

keep = pd_all$cell_id %in% cell_BC_assigned$cell_id
pd = pd_all[keep,]
count = count_all[,keep]
pd = pd[,c("cell_id","experiment_id")] %>% left_join(cell_BC_assigned[,c("cell_id","clonotype_id")], by = "cell_id")

### 51664 cells from 154 clonotypes

clonotype_list = unique(as.vector(pd$clonotype_id))
count_aggr = NULL
for(i in 1:length(clonotype_list)){
    print(paste0(i, "/", length(clonotype_list)))
    count_aggr = cbind(count_aggr, Matrix::rowSums(count[,pd$clonotype_id == clonotype_list[i]]))
}
colnames(count_aggr) = clonotype_list

exp = t(t(count_aggr) / colSums(count_aggr)) * 10000
exp = log(exp + 1)

gene_list = c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1")

for(i in 1:length(gene_list)){
    print(i)
    gene_i = gene_list[i]
    gene_id = mouse_gene$gene_ID[mouse_gene$gene_short_name == gene_i]
    
    dat = data.frame(exp = as.vector(exp[gene_id,]), clonotype_id = clonotype_list)
    dat = dat %>% left_join(hooke_ccs, by = "clonotype_id")
    target_regroup = if_else(dat$target == gene_i, gene_i, "NTC_and_Others")
    dat$target = factor(target_regroup, levels = c(gene_i, "NTC_and_Others"))
    
    fit = wilcox.test(dat$exp[dat$target == gene_i], dat$exp[dat$target == "NTC_and_Others"], alternative = "less")
    
    p = ggplot(dat, aes(target, exp, fill = target)) + geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2) + 
            labs(x="", y="Log normalized gene expression", title = paste0(gene_i, ", P-val = ", sprintf("%.2e", fit$p.val))) +
            theme_classic(base_size = 12) +
            theme(legend.position="none") +
            scale_fill_brewer(palette="Paired") +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/boxplot_gene_exp_", gene_i, "_NTC_and_Others.pdf"), p,
                   height  = 4, 
                   width = 3)
}





#############
### Pooled

hooke_ccs = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep2/hooke_ccs_pd.rds"))[[1]]
hooke_ccs = data.frame(clonotype_id = hooke_ccs$clonotype_id, target = hooke_ccs$target)
cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep2.rds"))
cell_BC_assigned = cell_BC_assigned[cell_BC_assigned$clonotype_id %in% hooke_ccs$clonotype_id,]

keep = pd_all$cell_id %in% cell_BC_assigned$cell_id
pd = pd_all[keep,]
count = count_all[,keep]
pd = pd[,c("cell_id","experiment_id")] %>% left_join(cell_BC_assigned[,c("cell_id","clonotype_id")], by = "cell_id")

### 13141 cells from 71 clonotypes

clonotype_list = unique(as.vector(pd$clonotype_id))
count_aggr = NULL
for(i in 1:length(clonotype_list)){
    print(paste0(i, "/", length(clonotype_list)))
    count_aggr = cbind(count_aggr, Matrix::rowSums(count[,pd$clonotype_id == clonotype_list[i]]))
}
colnames(count_aggr) = clonotype_list

exp = t(t(count_aggr) / colSums(count_aggr)) * 10000
exp = log(exp + 1)

gene_list = c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex")

for(i in 1:length(gene_list)){
    print(i)
    gene_i = gene_list[i]
    gene_id = mouse_gene$gene_ID[mouse_gene$gene_short_name == gene_i]
    
    dat = data.frame(exp = as.vector(exp[gene_id,]), clonotype_id = clonotype_list)
    dat = dat %>% left_join(hooke_ccs, by = "clonotype_id")
    target_regroup = if_else(dat$target == gene_i, gene_i, "NTC_and_Others")
    dat$target = factor(target_regroup, levels = c(gene_i, "NTC_and_Others"))
    
    fit = wilcox.test(dat$exp[dat$target == gene_i], dat$exp[dat$target == "NTC_and_Others"], alternative = "less")
    
    p = ggplot(dat, aes(target, exp, fill = target)) + geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2) + 
        labs(x="", y="Log normalized gene expression", title = paste0(gene_i, ", P-val = ", sprintf("%.2e", fit$p.val))) +
        theme_classic(base_size = 12) +
        theme(legend.position="none") +
        scale_fill_brewer(palette="Paired") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep2/boxplot_gene_exp_", gene_i, "_NTC_and_Others.pdf"), p,
           height  = 4, 
           width = 3)
}





