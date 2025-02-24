
######################################################################
### Identifying cell-type composition changes, on "arrayed" experiment
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

library(PLNmodels, lib.loc = "/net/gs/vol1/home/cxqiu/R/x86_64-pc-linux-gnu-library/4.3")
library(hooke)

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

#######################################################
### Step-2: Making separated UMAP on individual targets

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))

target_list = unique(cell_target$target)

for(i in target_list){
    print(i)
    clonotype_id_include = as.vector(unique(cell_target$clonotype_id[cell_target$target == i]))
    cell_include = cell_BC_assigned$cell_id[cell_BC_assigned$clonotype_id %in% clonotype_id_include]
    try(ggplot() +
            geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size=0.1) +
            geom_point(data = pd %>% filter(cell_id %in% cell_include), aes(x = UMAP_1, y = UMAP_2), size=0.5, color = "red") +
            theme_void() +
            theme(legend.position="none") +
            theme(plot.title = element_text(hjust = 0.5)) +
            ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/UMAP_rep1_",  i, ".png"),
                   dpi = 300,
                   height  = 5, 
                   width = 5), silent = TRUE)
}


df = cell_BC_assigned %>% 
    left_join(pd %>% select(cell_id, celltype), by = "cell_id") %>%
    left_join(cell_target, by = "clonotype_id") %>%
    group_by(target, celltype) %>%
    tally()
df$target = factor(df$target, levels = c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1","NONTARGETING"))
df$celltype = factor(df$celltype, levels = c("Primordial germ cells",
                                             "Epiblast",
                                             "Primitive streak",
                                             "Neuroectoderm",
                                             "Surface ectoderm",
                                             "Endoderm",
                                             "Paraxial mesoderm",
                                             "Lateral plate mesoderm",
                                             "Cardiomyocytes",
                                             "Endothelial cells",
                                             "Blood progenitors",
                                             "Erythroid cells"))

# Stacked + percent across clonotypes
p = df %>%
    ggplot(aes(fill=celltype, y=n, x=target)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=EB_celltype_color_code) +
    labs(x="", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/celltype_composition_target.pdf"),
           dpi = 300,
           height  = 5, 
           width = 5)


#####################################################################
### Step-3: performing hooker analysis to identify changed cell types

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
obj = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))

count = GetAssayData(obj, slot = "counts")[,as.vector(cell_BC_assigned$cell_id)]
pd_sub = pd[as.vector(cell_BC_assigned$cell_id),]
pd_sub$clonotype_id = as.vector(cell_BC_assigned$clonotype_id)
pd_sub = pd_sub %>% select(cell_id, experiment_id, UMAP_1, UMAP_2, celltype, clonotype_id) %>% 
    left_join(cell_target, by = "clonotype_id") %>% as.data.frame()
rownames(pd_sub) = as.vector(pd_sub$cell_id)

cds = new_cell_data_set(count,
                        cell_metadata = pd_sub)
cds = preprocess_cds(cds, use_genes = VariableFeatures(obj))
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_1", "UMAP_2")])

ccs = new_cell_count_set(cds, 
                         sample_group = "clonotype_id", 
                         cell_group = "celltype")

### Carm1        Gata6        Hand1         Hes5         Hhex         Lhx1
### 16           14            8           14           17           16
### Lmo2 NONTARGETING         Six3            T
### 11           19           18           21

ccm  = new_cell_count_model(ccs,
                            main_model_formula_str = "~target")

res = NULL
for (target_i in c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1")){
    print(target_i)
    cond_exp = estimate_abundances(ccm, tibble::tibble(target = target_i))
    cond_not_exp = estimate_abundances(ccm, tibble::tibble(target = "NONTARGETING"))
    
    cond_ne_v_e_tbl = compare_abundances(ccm, cond_not_exp, cond_exp)

    res = rbind(res, cond_ne_v_e_tbl %>% select(cell_group, target_x, target_y,
                                                delta_log_abund, delta_log_abund_se, delta_p_value, delta_q_value))
}
res$fdr = p.adjust(res$delta_p_value, method = "fdr")
res %>% filter(fdr < 0.1) %>% arrange(target_y, delta_log_abund) %>% as.data.frame()
saveRDS(res, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/hooke_result.rds"))

write.table(res, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/hooke_result.txt"), row.names=F, quote=F, sep="\t")

##################################################
### Step-4:making heatmap based on delta_log_abund

dat = res %>% select(target_y, cell_group, delta_log_abund) %>%
    dcast(target_y ~ cell_group)
rownames(dat) = dat[,1]
dat = dat[,-1]

dat = dat[c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1"),
          c("Primordial germ cells",
            "Epiblast",
            "Primitive streak",
            "Neuroectoderm",
            "Surface ectoderm",
            "Endoderm",
            "Paraxial mesoderm",
            "Lateral plate mesoderm",
            "Cardiomyocytes",
            "Endothelial cells",
            "Blood progenitors",
            "Erythroid cells")]

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/hooke_heatmap.pdf"), 8, 5)
heatmap.2(as.matrix(t(dat)), 
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

dat = res %>% filter(fdr < 0.1) %>% select(target_y, cell_group, delta_log_abund) %>%
    dcast(target_y ~ cell_group)
rownames(dat) = dat[,1]
dat = dat[,-1]

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/hooke_heatmap_only_sig.pdf"), 8, 5)
heatmap.2(as.matrix(t(dat)), 
          col=Colors, 
          scale="none", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()


############################################################
### Step-5: making boxplot with cell-type fractions changing

res_sig = res %>% filter(fdr < 0.1)

cell_freq = exprs(ccs)
cell_freq = 100*t(t(cell_freq)/colSums(cell_freq))

for(i in 1:nrow(res_sig)){
    target_i = res_sig$target_y[i]
    celltype_i = res_sig$cell_group[i]
    cell_freq_target = cell_freq[rownames(cell_freq) == celltype_i,
                                 ccs$target == target_i]
    cell_freq_control = cell_freq[rownames(cell_freq) == celltype_i,
                                 ccs$target == "NONTARGETING"]
    df = data.frame(frac = c(as.vector(cell_freq_target), as.vector(cell_freq_control)),
                    group = c(rep(target_i, times = length(cell_freq_target)), rep("NONTARGETING", times = length(cell_freq_control))))
    
    df$group = factor(df$group, levels = c(target_i, "NONTARGETING"))
    
    try(ggplot(df, aes(group, frac, fill = group)) + geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2) + 
        labs(x="", y="% of cells", title = celltype_i) +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        scale_fill_brewer(palette="Paired") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/boxplot_", target_i, "_", gsub(" ","_",celltype_i), ".pdf"),
               height  = 5, 
               width = 3), silent = T)
    
}


### performing wilcox-test for Carm1
x_1 = as.vector(cell_freq[rownames(cell_freq) == "Epiblast",
                          ccs$target == "Carm1"])
x_2 = as.vector(cell_freq[rownames(cell_freq) == "Epiblast",
                          ccs$target == "NONTARGETING"])
wilcox.test(x_1, x_2) ### p = 0.01

x_1 = as.vector(cell_freq[rownames(cell_freq) == "Primordial germ cells",
                          ccs$target == "Carm1"])
x_2 = as.vector(cell_freq[rownames(cell_freq) == "Primordial germ cells",
                          ccs$target == "NONTARGETING"])
wilcox.test(x_1, x_2) ### p = 0.084


x_1 = as.vector(cell_freq[rownames(cell_freq) == "Cardiomyocytes",
                          ccs$target == "Hand1"])
x_2 = as.vector(cell_freq[rownames(cell_freq) == "Cardiomyocytes",
                          ccs$target == "NONTARGETING"])
wilcox.test(x_1, x_2) ### p = 0.0001694

x_1 = as.vector(cell_freq[rownames(cell_freq) == "Lateral plate mesoderm",
                          ccs$target == "Hand1"])
x_2 = as.vector(cell_freq[rownames(cell_freq) == "Lateral plate mesoderm",
                          ccs$target == "NONTARGETING"])
wilcox.test(x_1, x_2) ### p = 0.0009234


###########################################
### Step-6: plot gene expression of targets

cds = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cds.rds"))

gene_list = c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1")

for(gene_i in gene_list){
    print(gene_i)
    try(my_plot_cells(cds, genes = gene_i, cell_size = 0.5) + 
        theme_void() + NoLegend() + theme(strip.text.x = element_blank()) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/gene_expression_", gene_i, ".png"), 
               dpi = 300, height = 5, width = 5), silent = T)
}


####################################################################################
### Step-7: making a big heatmap, each column is a clonotype, each row is a celltype 

cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep1.rds"))
cell_BC_assigned_num = cell_BC_assigned %>% group_by(clonotype_id) %>% tally() %>% filter(n >= 100)
cell_BC_assigned = cell_BC_assigned[cell_BC_assigned$clonotype_id %in% as.vector(cell_BC_assigned_num$clonotype_id),]

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))

pd_sub = pd[as.vector(cell_BC_assigned$cell_id),]
pd_sub$clonotype_id = as.vector(cell_BC_assigned$clonotype_id)
pd_sub = pd_sub %>% select(cell_id, experiment_id, UMAP_1, UMAP_2, celltype, clonotype_id) %>% 
    left_join(cell_target, by = "clonotype_id") %>% as.data.frame()
rownames(pd_sub) = as.vector(pd_sub$cell_id)

dat = pd_sub %>% group_by(clonotype_id, celltype) %>% tally() %>%
    dcast(celltype~clonotype_id, fill = 0)
rownames(dat) = as.vector(dat[,1])
dat = as.matrix(dat[,-1])

dat = dat/apply(dat, 1, sum)
dat_norm = t(dat)/apply(dat, 2, sum)

col_names = c("Primordial germ cells",
  "Epiblast",
  "Primitive streak",
  "Neuroectoderm",
  "Surface ectoderm",
  "Endoderm",
  "Paraxial mesoderm",
  "Lateral plate mesoderm",
  "Cardiomyocytes",
  "Endothelial cells",
  "Blood progenitors",
  "Erythroid cells")

row_names = NULL
for(i in col_names){
    for(j in 1:nrow(dat_norm)){
        tmp = names(sort(dat_norm[j,], decreasing=T))[1]
        if(tmp == i & !rownames(dat_norm)[j] %in% row_names){
            row_names = c(row_names, rownames(dat_norm)[j])
        }
    }
}
dat_norm = dat_norm[row_names,col_names]

clonotypeID = rownames(dat_norm)
cell_target_x = data.frame(clonotype_id = clonotypeID) %>%
    left_join(cell_target, by = "clonotype_id")

library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/clonotype_celltype_composition_heatmap_rep1.pdf"), 8, 5)
h = heatmap.2(t(dat_norm), 
          col=Colors, 
          scale="col", 
          Rowv = T, 
          Colv = T, 
          ColSideColors=TF_screen_color_code[as.vector(cell_target_x$target)],
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 1, 
          cexCol = 1,
          margins = c(5,5))
dev.off()

clonotypeID = rownames(dat_norm)[h$colInd]
cell_target_y = data.frame(clonotype_id = clonotypeID) %>%
    left_join(cell_target, by = "clonotype_id")



###############################################
### Step-8: making UMAP on individual clonotype

cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep1.rds"))
pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
pd_sub = pd[as.vector(cell_BC_assigned$cell_id),]
pd_sub$clonotype_id = as.vector(cell_BC_assigned$clonotype_id)

target_list = unique(cell_target$target)
sort(table(cell_target$target))

for(target_i in target_list){
    print(target_i)
    
    n = sum(cell_target$target == target_i)
    col_num = 10
    if(n %% col_num == 0){
        row_num = n/col_num
    } else (
        row_num = floor(n/col_num) + 1
    )
    
    try(ggplot() +
            geom_point(data = pd_sub %>%
                           as.data.frame() %>%
                           dplyr::select(-clonotype_id),
                       aes(x = UMAP_1,
                           y = UMAP_2),
                       color = "grey80",
                       stroke = 0,
                       size = 0.6) +
            geom_point(data = pd_sub %>%
                           as.data.frame() %>%
                           filter(clonotype_id %in% as.vector(cell_target$clonotype_id[cell_target$target == target_i])),
                       aes(x = UMAP_1,
                           y = UMAP_2),
                       color = TF_screen_color_code[target_i],
                       stroke = 0,
                       size = 0.6) +
            monocle3:::monocle_theme_opts() +
            facet_wrap(~clonotype_id, ncol = col_num, nrow = row_num) +
            theme_void() +
            theme(legend.position = "none",
                  strip.text = element_blank()) + 
            ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/UMAP_clonotype_rep1_", target_i, ".png"),
                   dpi = 300,
                   height  = row_num, 
                   width = 10), silent = T)
    
}






###############################
### Step-9: pseudo-bulk samples

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
obj = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))

count = GetAssayData(obj, slot = "counts")[,as.vector(cell_BC_assigned$cell_id)]
pd_sub = pd[as.vector(cell_BC_assigned$cell_id),]
pd_sub$clonotype_id = as.vector(cell_BC_assigned$clonotype_id)
pd_sub = pd_sub %>% select(cell_id, experiment_id, UMAP_1, UMAP_2, celltype, clonotype_id) %>% 
    left_join(cell_target, by = "clonotype_id") %>% as.data.frame()
rownames(pd_sub) = as.vector(pd_sub$cell_id)

clonotype_list = as.vector(cell_target$clonotype_id)
count_aggr = NULL
for(i in 1:length(clonotype_list)){
    print(paste0(i, "/", length(clonotype_list)))
    count_aggr = cbind(count_aggr, Matrix::rowSums(count[,pd_sub$clonotype_id == clonotype_list[i]]))
}


### pseudo-bulk PCA (START) ###
obj_aggr = CreateSeuratObject(count_aggr, meta.data = cell_target)
obj_aggr = NormalizeData(obj_aggr, normalization.method = "LogNormalize", scale.factor = 10000)
obj_aggr = FindVariableFeatures(obj_aggr, selection.method = "vst", nfeatures = 3000)
gene_use = VariableFeatures(obj_aggr)

cds_aggr = new_cell_data_set(count_aggr,
                             cell_metadata = cell_target)

set.seed(2016)
FM = monocle3:::normalize_expr_data(cds_aggr, 
                                    norm_method = "log", 
                                    pseudo_count = 1)
FM = FM[gene_use,]

num_dim = 10
scaling = TRUE
set.seed(2016)
irlba_res = my_sparse_prcomp_irlba(Matrix::t(FM), 
                                   n = min(num_dim, min(dim(FM)) - 1), 
                                   center = scaling, 
                                   scale. = scaling)
preproc_res = irlba_res$x

prop_var_expl = irlba_res$sdev^2/sum(irlba_res$sdev^2)
print(prop_var_expl)
# 0.28294665 0.19023439 0.11288998 0.08620485 0.07977272 0.06330701
# 0.05737338 0.04557207 0.04292819 0.03877076

df = data.frame(clonotype_id = as.vector(cds_aggr$clonotype_id),
                PC_1 = preproc_res[,1],
                PC_2 = preproc_res[,2],
                PC_3 = preproc_res[,3],
                target = as.vector(cds_aggr$target))

save_path = "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup"
fig = plot_ly(df, x=~PC_1, y=~PC_2, z=~PC_3, size = I(300), color = ~target, colors = TF_screen_color_code)
saveWidget(fig, paste0(save_path, "/sam_tf/Hand_picked_pseudobulk_PCA.html"), selfcontained = FALSE, libdir = "tmp")
### pseudo-bulk PCA (DONE) ###


### gene expression across cells with NTC and target gRNAs
### normalization on each column
exp = t(t(count) / colSums(count)) * 10000
exp@x = log(exp@x + 1)

gene_list = c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1")

df = NULL

for(i in 1:length(gene_list)){
    print(i)
    gene_i = gene_list[i]
    gene_id = mouse_gene$gene_ID[mouse_gene$gene_short_name == gene_i]
    
    df_i = data.frame(exp = as.vector(exp[gene_id,]),
                    target = as.vector(pd_sub$target)) 
    
    df = rbind(df, data.frame(gene = gene_i,
               mean_exp_ko = mean(df_i$exp[df_i$target == gene_i]),
               mean_exp_ntc = mean(df_i$exp[df_i$target == "NONTARGETING"])))
}

df %>% 
    ggplot(aes(mean_exp_ntc, mean_exp_ko, color = gene)) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_text(aes(label = gene), color = "black") +
    labs(x="Mean log-normalized gene expression in cells with NTC", y="Mean log-normalized gene expression in cells with target gRNAs", title="") +
    theme_classic(base_size = 10) +
    scale_color_manual(values=TF_screen_color_code) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/gene_exp_compare.pdf"), 
           height  = 5, 
           width = 5)



### boxplot between individual EBs

exp = t(t(count_aggr) / colSums(count_aggr)) * 10000
exp = log(exp + 1)

gene_list = c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1")

for(i in 1:length(gene_list)){
    print(i)
    gene_i = gene_list[i]
    gene_id = mouse_gene$gene_ID[mouse_gene$gene_short_name == gene_i]
    
    df = data.frame(exp = as.vector(exp[gene_id,]),
                    target = as.vector(cell_target$target)) %>% filter(target %in% c(gene_i, "NONTARGETING"))
    df$target = factor(df$target, levels = c(gene_i, "NONTARGETING"))
    
    try(ggplot(df, aes(target, exp, fill = target)) + geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2) + 
            labs(x="", y="Log-normalized gene expression", title = gene_i) +
            theme_classic(base_size = 10) +
            theme(legend.position="none") +
            scale_fill_brewer(palette="Paired") +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
            ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/boxplot_gene_exp_", gene_i, ".pdf"),
                   height  = 5, 
                   width = 3), silent = T)
}

### boxplot between individual EBs, separated by celltypes

cell_target_celltype = pd_sub %>% group_by(clonotype_id, target, celltype) %>% tally()
count_aggr_celltype = NULL
for(i in 1:nrow(cell_target_celltype)){
    print(paste0(i, "/", nrow(cell_target_celltype)))
    count_aggr_celltype = cbind(count_aggr_celltype, 
                                Matrix::rowSums(count[,pd_sub$clonotype_id == cell_target_celltype$clonotype_id[i] &
                                                          pd_sub$celltype == cell_target_celltype$celltype[i], drop=FALSE]))
}


exp = t(t(count_aggr_celltype) / colSums(count_aggr_celltype)) * 10000
exp = log(exp + 1)

gene_list = c("Hand1")

for(i in 1:length(gene_list)){
    print(i)
    gene_i = gene_list[i]
    gene_id = mouse_gene$gene_ID[mouse_gene$gene_short_name == gene_i]
    
    df = data.frame(exp = as.vector(exp[gene_id,]),
                    target = as.vector(cell_target_celltype$target),
                    celltype = as.vector(cell_target_celltype$celltype)) %>% filter(target %in% c(gene_i, "NONTARGETING"))
    df$target = factor(df$target, levels = c(gene_i, "NONTARGETING"))
    
    try(ggplot(df, aes(target, exp, fill = target)) + geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2) + 
            ggh4x::facet_grid2(.~celltype, scales = "free_y", independent = "y") +
            labs(x="", y="Log-normalized gene expression", title = gene_i) +
            theme_classic(base_size = 8) +
            theme(legend.position="none") +
            scale_fill_brewer(palette="Paired") +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.text.x = element_text(color="black", angle = 90), axis.text.y = element_text(color="black")) +
            ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/boxplot_gene_exp_by_celltype_", gene_i, ".pdf"),
                   height  = 3, 
                   width = 12), silent = T)
}



###############################################################################
### Step-10: within the neighbors, % of the same gRNAs vs. mean gene expression

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
obj = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))

count = GetAssayData(obj, slot = "counts")[,as.vector(cell_BC_assigned$cell_id)]
pca_coor = Embeddings(obj, reduction = "pca")[as.vector(cell_BC_assigned$cell_id),]
pd_sub = pd[as.vector(cell_BC_assigned$cell_id),]
pd_sub$clonotype_id = as.vector(cell_BC_assigned$clonotype_id)
pd_sub = pd_sub %>% select(cell_id, experiment_id, UMAP_1, UMAP_2, celltype, clonotype_id) %>% 
    left_join(cell_target, by = "clonotype_id") %>% as.data.frame()
rownames(pd_sub) = as.vector(pd_sub$cell_id)

kadj = 20
k.param = kadj + 1; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = pca_coor,
    k = k.param,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)
nn.ranked = Indices(object = nn.ranked)
nn_matrix = nn.ranked[,-1]
nn_matrix = as.matrix(nn_matrix)

dat_gRNA = NULL
for(kk in 1:kadj){
    dat_gRNA = cbind(dat_gRNA, pd_sub$target[nn_matrix[,kk]])
}

exp = t(t(count) / colSums(count)) * 10000
exp@x = log(exp@x + 1)

gene_list = c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1")

for(i in 1:length(gene_list)){
    gene_i = gene_list[i]; print(gene_i)
    gene_id = mouse_gene$gene_ID[mouse_gene$gene_short_name == gene_i]
    
    pd_sub$gene_exp = as.vector(exp[gene_id,])
    
    dat_exp = NULL
    for(kk in 1:kadj){
        dat_exp = cbind(dat_exp, pd_sub$gene_exp[nn_matrix[,kk]])
    }
    
    df = data.frame(mean_freq_gRNA = apply(dat_gRNA == gene_i, 1, mean),
                    mean_gene_exp = apply(dat_exp, 1, mean))
    
    fit = cor.test(df$mean_freq_gRNA, df$mean_gene_exp)
    
    df$mean_freq_gRNA = factor(df$mean_freq_gRNA)
    
    try(df %>% 
        ggplot(aes(mean_freq_gRNA, mean_gene_exp)) + geom_boxplot() +
        labs(x="% of neighbors from target TF", y="Mean log-normalized gene expression of neighbors", title=paste0(gene_i, ", Pearson corr = ", round(fit$estimate, 3), ", pval = ", round(fit$p.value, 3))) +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/Neighbors_", gene_i, ".pdf"), 
               dpi = 300,
               height  = 5, 
               width = 5), silent = T)
}








#############################################################################################
### Step-11: scatterplot of TF expression vs cell type proportion per EB (for each cell type)

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
obj = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))

count = GetAssayData(obj, slot = "counts")[,as.vector(cell_BC_assigned$cell_id)]
pca_coor = Embeddings(obj, reduction = "pca")[as.vector(cell_BC_assigned$cell_id),]
pd_sub = pd[as.vector(cell_BC_assigned$cell_id),]
pd_sub$clonotype_id = as.vector(cell_BC_assigned$clonotype_id)
pd_sub = pd_sub %>% select(cell_id, experiment_id, UMAP_1, UMAP_2, celltype, clonotype_id) %>% 
    left_join(cell_target, by = "clonotype_id") %>% as.data.frame()
rownames(pd_sub) = as.vector(pd_sub$cell_id)

exp = t(t(count) / colSums(count)) * 10000
exp@x = log(exp@x + 1)

gene_list = c("Gata6")

for(i in 1:length(gene_list)){
    gene_i = gene_list[i]; print(gene_i)
    gene_id = mouse_gene$gene_ID[mouse_gene$gene_short_name == gene_i]
    
    pd_sub$gene_exp = as.vector(exp[gene_id,])
    
    df = pd_sub %>% filter(target == gene_i) %>%
        group_by(clonotype_id, celltype) %>%
        summarize(mean_exp = mean(gene_exp), cell_num = n()) %>%
        left_join(pd_sub %>% filter(target == gene_i) %>%
                      group_by(clonotype_id) %>% tally() %>% rename(cell_num_total = n), by = "clonotype_id") %>%
        mutate(frac = 100*cell_num/cell_num_total) %>% select(clonotype_id, celltype, mean_exp, frac)
    
    try(ggplot(df, aes(mean_exp, frac, color = clonotype_id)) + geom_point() +
            ggh4x::facet_grid2(.~celltype, scales = "free", independent = "all") +
            labs(x="Mean log-normalized gene expression", y="% of cells", title = gene_i) +
            theme_classic(base_size = 8) +
            theme(legend.position="none") +
            scale_fill_brewer(palette="Paired") +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.text.x = element_text(color="black", angle = 90), axis.text.y = element_text(color="black")) +
            ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/TF_exp_vs_celltype_frac_", gene_i, ".pdf"),
                   height  = 3, 
                   width = 20), silent = T)
    
}







