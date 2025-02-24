
######################################################
### Comparing cell-type-compositions between KO and WT
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tf
### big_screen.cell_gene.rds

pd = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned_pd.rds"))
pd$cell = gsub("all_TF_screen_GEx_", "", as.vector(pd$cell_id))
cell_gene = readRDS("big_screen.cell_gene.rds")

######################################################################################################
### Step-1: regroup cells, to ensure more cells in each cell cluster (only for testing, not necessary)

celltype_regroup = rep(NA, nrow(pd))

celltype_regroup[pd$celltype %in% c("Primordial germ cells",
                                    "ESCs (2-cell state)",
                                    "Epiblast",
                                    "Epiblast and primitive streak")] = "Epiblast"
celltype_regroup[pd$celltype %in% c("Primitive streak")] = "Primitive streak"
celltype_regroup[pd$celltype %in% c("Neuroectoderm",
                                    "Early neurons",
                                    "Floor plate",
                                    "Eye field")] = "Neuroectoderm"
celltype_regroup[pd$celltype %in% c("Surface ectoderm")] = "Surface ectoderm"
celltype_regroup[pd$celltype %in% c("Notochord",
                                    "Definitive endoderm")] = "Notochord and def. endoderm"
celltype_regroup[pd$celltype %in% c("Extraembryonic visceral endoderm")] = "ExE visceral endoderm"
celltype_regroup[pd$celltype %in% c("Parietal endoderm")] = "Parietal endoderm"
celltype_regroup[pd$celltype %in% c("Nascent mesoderm",
                                    "Lateral plate mesoderm",
                                    "Cardiomyocytes")] = "Mesoderm"
celltype_regroup[pd$celltype %in% c("Endothelial cells",
                                    "Blood progenitors",
                                    "Erythroid cells")] = "Blood"

pd$celltype_regroup = as.vector(celltype_regroup) 


cell_gene = cell_gene %>% left_join(pd[,c("cell","celltype","celltype_regroup")], by = "cell")


##########################################################
### Step-2: filtering out targets with low number of cells

### this is the target level analysis, ignoring the heterogenity between guides within the same target

target_cell_num = cell_gene %>% filter(gene != "NONTARGETING") %>% 
    select(cell, gene) %>% unique() %>%
    group_by(gene) %>% tally()

NTC_guide_num = cell_gene %>% filter(gene == "NONTARGETING") %>% 
    select(cell, barcode) %>% unique() %>%
    group_by(barcode) %>% tally()

### filtering out targets with less than 50 cells

target_cell_num = target_cell_num %>% filter(n >= 50)
print(nrow(target_cell_num)) ### n = 820

### filtering out NTCs with less than 20 cells

NTC_guide_num = NTC_guide_num %>% filter(n >= 20)
print(nrow(NTC_guide_num)) ### n = 540




##############################################################
### Step-3: performing chi-squared test on each target vs. ntc

mat = cell_gene %>% group_by(gene, celltype) %>% 
    tally() %>% dcast(gene~celltype, fill = 0)
rownames(mat) = as.vector(mat$gene)
mat = mat[,-1]
mat = as.matrix(mat)

res = NULL
for(i in 1:nrow(target_cell_num)){
    target_i = as.vector(target_cell_num$gene)[i]
    contingency_table = matrix(c(mat[target_i, ], mat["NONTARGETING", ]), nrow = 2, byrow = TRUE)
    chi_sq_test = chisq.test(contingency_table)
    res = c(res, chi_sq_test$p.value)
}

res_x = data.frame(gene = as.vector(target_cell_num$gene),
                 cell_num = as.vector(target_cell_num$n),
                 pval = res)

res = data.frame(gene = as.vector(target_cell_num$gene),
                 pval = res)


### quick checking if the siginificant p-values are correlated with the cell number for individual targets
res_x$log10pval = -log10(res_x$pval)
res_x$log10pval[is.infinite((res_x$log10pval))] = 320
res_x$log2_cell_num = log2(res_x$cell_num)

cor.test(res_x$log2_cell_num, res_x$log10pval, method = "spearman")
cor.test(res_x$log2_cell_num[res_x$log10pval > 100], res_x$log10pval[res_x$log10pval > 100], method = "spearman")

p = ggplot() + 
    geom_point(data = res_x, aes(x = log10pval, y = log2_cell_num)) + 
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/cell_num_chi_squared_pval.pdf"))
print(p)
dev.off()



######################################################################
### Step-4: performing chi-squared test on each two NTC guides vs. ntc


mat = cell_gene %>% filter(gene == "NONTARGETING") %>% group_by(barcode, celltype) %>% 
    tally() %>% dcast(barcode~celltype, fill = 0)
rownames(mat) = as.vector(mat$barcode)
mat = mat[,-1]
mat = as.matrix(mat)

mat_sum = apply(mat, 2, sum)

res_NTC = NULL
for(iter in 1:1){
    print(iter)
    
    each_group = 3
    group_num = floor(nrow(NTC_guide_num)/each_group)
    set.seed(1234)
    NTC_guide_num$group = sample(rep(1:group_num, each = each_group))[1:nrow(NTC_guide_num)]
    
    for(i in 1:group_num){
        barcode_include = as.vector(NTC_guide_num$barcode[NTC_guide_num$group == i])
        
        mat_sub = cell_gene %>% filter(barcode %in% barcode_include) %>% 
            slice_sample(n = 100) %>%
            group_by(barcode, celltype) %>% 
            tally() %>% dcast(barcode~celltype, fill = 0)
        rownames(mat_sub) = as.vector(mat_sub$barcode)
        mat_sub = apply(as.matrix(mat_sub[,-1]), 2, sum)
        mat_sub = mat_sub[names(mat_sum)]
        mat_sub[is.na(mat_sub)] = 0
        
        contingency_table = matrix(c(mat_sub, mat_sum), nrow = 2, byrow = TRUE)
        chi_sq_test = chisq.test(contingency_table)
        res_NTC = c(res_NTC, chi_sq_test$p.value)
    }
}

res_sig = res[res$pval < quantile(res_NTC, 0.05),]
res_sig = res_sig[order(res_sig$pval),]
saveRDS(list(res, res_NTC, res_sig), paste0(work_path, "/analysis/all_TF_screen_GEx/chi_squared_test.rds"))


############################
### Step-5: making a QQ plot

res_combine = data.frame(pval = c(res$pval, res_NTC),
                         log10pval = c(-log10(res$pval), -log10(res_NTC)),
                         group = c(rep("TF", nrow(res)), rep("NTC", length(res_NTC))))
res_combine$log10pval[is.infinite((res_combine$log10pval))] = 320

p = res_combine %>% ggplot(aes(sample = log10pval, color = group)) + 
    stat_qq(distribution = qexp, dparams=list(rate=log(10)),size = 2) +
    geom_hline(yintercept = -log10(quantile(res_NTC, 0.05)), linetype="dotted") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/chi_squared_test_QQplot_x.pdf"))
print(p)
dev.off()

qlog10 = function(dat){
    # assuming "nn_de" is a dataframe where each row contains one of your p-values
    n <- length(dat)
    
    # Generate quantiles from a normal distribution
    quantiles <- qnorm(seq(0, 1, length.out = n + 1)[-1])
    
    # Rescale quantiles to fit between 0.00001 and 1
    perfect_values <- pmin(pmax(pnorm(quantiles), 1e-5), 1)
    
    # Apply -log10() function to each value
    log_dis <- -log10(perfect_values)
    
    expected_p <- log_dis[order(-log_dis)]
    
    return(expected_p)
}

dat_1 = subset(res_combine, group == "TF")
dat_1 = dat_1[order(dat_1$pval),]
dat_1$expected_p = qlog10(dat_1$pval)

dat_2 = subset(res_combine, group == "NTC")
dat_2 = dat_2[order(dat_2$pval),]
dat_2$expected_p = qlog10(dat_2$pval)

p = ggplot() + 
    geom_point(data = dat_1, aes(x = expected_p, y = log10pval), color = "blue") + 
    geom_point(data = dat_2, aes(x = expected_p, y = log10pval), color = "red") + 
    geom_hline(yintercept = -log10(quantile(res_NTC, 0.05)), linetype="dotted") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/chi_squared_test_QQplot.pdf"))
print(p)
dev.off()


####################################################################################
### Step-6: calculating odds ratio for the top 3 TF targets for individual cell type

mat = cell_gene %>% group_by(gene, celltype) %>% 
    tally() %>% dcast(gene~celltype, fill = 0)
rownames(mat) = as.vector(mat$gene)
mat = mat[,-1]
mat = as.matrix(mat)

dat_res = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/chi_squared_test.rds"))
res_sig = dat_res[[1]]  ### n = 820 TF targets

mat = mat[rownames(mat) %in% c(as.vector(res_sig$gene), "NONTARGETING"),]

### top significant TF targets
res_or = NULL
res_pval = NULL
top_n = nrow(res_sig)
for(i in 1:top_n){
    print(paste0(i, "/", top_n))
    or_i = NULL
    pval_i = NULL
    mat_1 = mat[res_sig$gene[i],]; mat_1_sum = sum(mat_1)
    mat_2 = mat["NONTARGETING",]; mat_2_sum = sum(mat_2)
    for(j in 1:ncol(mat)){
        a = mat_1[j]; b = mat_1_sum - a
        c = mat_2[j]; d = mat_2_sum - c
        x = matrix(c(a,b,c,d), 2, 2)
        fit = fisher.test(x)
        or_i = c(or_i, fit$estimate)
        pval_i = c(pval_i, fit$p.value)
    }
    res_or = rbind(res_or, or_i)
    res_pval = rbind(res_pval, pval_i)
}
rownames(res_or) = rownames(res_pval) = as.vector(res_sig$gene[1:top_n])
colnames(res_or) = colnames(res_pval) = colnames(mat)


### plot heatmap for individual cell types with top three targets (order by odds ratio)

res = melt(as.matrix(res_or)) %>% as.data.frame() %>% rename(or = value) %>%
    left_join(melt(as.matrix(res_pval)) %>% as.data.frame() %>% rename(pval = value), by = c("Var1", "Var2")) %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>%
    filter(fdr < 0.05) %>%
    group_by(Var2) %>% slice_max(order_by = or, n = 3, with_ties = F)

aggr_matrix = res_or[rownames(res_or) %in% as.vector(res$Var1),]
# manually order rownames (cell types)
col_names = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% colnames(aggr_matrix)]
row_names = NULL
for(i in col_names){
    for(j in 1:nrow(aggr_matrix)){
        tmp = names(sort(aggr_matrix[j,], decreasing=T))[1]
        if(tmp == i & !rownames(aggr_matrix)[j] %in% row_names){
            row_names = c(row_names, rownames(aggr_matrix)[j])
        }
    }
}
aggr_matrix = aggr_matrix[row_names,col_names]



library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/top_3TFs_celltype.pdf"), 8, 5)
heatmap.2(as.matrix(t(aggr_matrix)), 
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

write.table(rownames(aggr_matrix), paste0(work_path, "/analysis/all_TF_screen_GEx/plot/top_3TFs_celltype.rownames.txt"), row.names=F, col.names=F, quote=F, sep="\t")


### plot heatmap for individual cell types with top three targets (order by p-value)

### filter target with at least 50 cells
target_cell_num = cell_gene %>% filter(gene != "NONTARGETING") %>% 
    select(cell, gene) %>% unique() %>%
    group_by(gene) %>% tally()

target_cell_num = target_cell_num %>% filter(n >= 50)
print(nrow(target_cell_num)) ### n = 820

res = melt(as.matrix(res_or)) %>% as.data.frame() %>% rename(or = value) %>%
    left_join(melt(as.matrix(res_pval)) %>% as.data.frame() %>% rename(pval = value), by = c("Var1", "Var2")) %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>%
    filter(Var1 %in% as.vector(target_cell_num$gene)) %>% 
    group_by(Var2) %>% slice_min(order_by = fdr, n = 3, with_ties = F)

aggr_matrix = res_or[rownames(res_or) %in% as.vector(res$Var1),]
# manually order rownames (cell types)
col_names = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% colnames(aggr_matrix)]
row_names = NULL
for(i in col_names){
    for(j in 1:nrow(aggr_matrix)){
        tmp = names(sort(aggr_matrix[j,], decreasing=T))[1]
        if(tmp == i & !rownames(aggr_matrix)[j] %in% row_names){
            row_names = c(row_names, rownames(aggr_matrix)[j])
        }
    }
}
aggr_matrix = aggr_matrix[row_names,col_names]



library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/top_3TFs_celltype_order_by_pval.pdf"), 8, 5)
heatmap.2(as.matrix(t(log2(aggr_matrix+0.01))), 
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

write.table(rownames(aggr_matrix), paste0(work_path, "/analysis/all_TF_screen_GEx/plot/top_3TFs_celltype_order_by_pval.rownames.txt"), row.names=F, col.names=F, quote=F, sep="\t")




###########################################
### Step-7: plot gRNA for individual target


pd = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned_pd.rds"))
cell_gene = readRDS(paste0(work_path, "/processing/all_TF_screen_sgRNA/cell_gene.rds"))
pd$cell = gsub("all_TF_screen_GEx_", "", as.vector(pd$cell_id))

cnt = 1
for(target_i in as.vector(target_cell_num$gene)){
    print(cnt)
    cnt = cnt + 1
    
    df = pd %>% select(cell, UMAP_1, UMAP_2) %>% 
        mutate(target = if_else(cell %in% as.vector(cell_gene$cell[cell_gene$gene == target_i]), target_i, "none"))
    
    try(ggplot() +
            geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
            geom_point(data = subset(df, target == target_i), aes(x = UMAP_1, y = UMAP_2), size=1, color = "red") +
            theme_void() +
            theme(legend.position="none") +
            labs(title = target_i) +
            theme(plot.title = element_text(hjust = 0.5)) +
            ggsave(paste0(work_path, "/analysis/all_TF_screen_GEx/sgRNA_plot/", target_i, ".png"),
                   dpi = 300,               
                   height  = 5, 
                   width = 5), silent = TRUE)
}





##################################################################################
### Step-8: investigating the heterogenity of different guides for the same target

target_cell_num = cell_gene %>% filter(gene != "NONTARGETING") %>% 
    select(cell, gene) %>% unique() %>%
    group_by(gene) %>% tally() %>% rename(target_num = n) %>% filter(target_num >= 50)

guide_cell_num = cell_gene %>% filter(gene != "NONTARGETING") %>% 
    select(cell, barcode, gene) %>% unique() %>%
    group_by(barcode, gene) %>% tally() %>% rename(guide_num = n) %>% 
    filter(gene %in% target_cell_num$gene) %>%
    left_join(target_cell_num, by = "gene") %>%
    mutate(frac = 100*guide_num/target_num)

guide_cell_num_top1 = guide_cell_num %>% group_by(gene) %>% slice_max(order_by = frac, n = 1, with_ties = F)
guide_cell_num_top2 = guide_cell_num %>% filter(!barcode %in% guide_cell_num_top1$barcode) %>% 
    group_by(gene) %>% slice_max(order_by = frac, n = 1, with_ties = F)
guide_cell_num_top3 = guide_cell_num %>% filter(!barcode %in% c(guide_cell_num_top1$barcode, guide_cell_num_top2$barcode)) %>% 
    group_by(gene) %>% slice_max(order_by = frac, n = 1, with_ties = F)

df = data.frame(frac = c(guide_cell_num_top1$frac, guide_cell_num_top2$frac, guide_cell_num_top3$frac),
                group = c(rep("top1", nrow(guide_cell_num_top1)), rep("top2", nrow(guide_cell_num_top2)), rep("top3", nrow(guide_cell_num_top3))))

df$group = factor(df$group, levels = paste0("top", c(3,2,1)))

p = ggplot(df, aes(group, frac, fill = group)) + geom_boxplot() +
    labs(x="", y="") +
    theme_classic(base_size = 10) +
    coord_flip() +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set2") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf("~/share/guide_frac.pdf",3,5)
print(p)
dev.off()


##############

guide_cell_num = cell_gene %>% filter(gene != "NONTARGETING") %>% 
    select(cell, barcode, gene) %>% unique() %>%
    group_by(barcode, gene) %>% tally() %>% filter(n >= 50)

mat = cell_gene %>% group_by(gene, celltype) %>% 
    tally() %>% dcast(gene~celltype, fill = 0)
rownames(mat) = as.vector(mat$gene)
mat = mat[,-1]
mat_gene = as.matrix(mat)

mat = cell_gene %>% group_by(barcode, celltype) %>% 
    tally() %>% dcast(barcode~celltype, fill = 0)
rownames(mat) = as.vector(mat$barcode)
mat = mat[,-1]
mat = as.matrix(mat)

res = NULL
for(i in 1:nrow(guide_cell_num)){
    guide_i = as.vector(guide_cell_num$barcode)[i]
    contingency_table = matrix(c(mat[guide_i, ], mat_gene["NONTARGETING", ]), nrow = 2, byrow = TRUE)
    chi_sq_test = chisq.test(contingency_table)
    res = c(res, chi_sq_test$p.value)
}
res = data.frame(barcode = as.vector(guide_cell_num$barcode),
                 gene = as.vector(guide_cell_num$gene),
                 pval = res)




guide_NTC = cell_gene %>% filter(gene == "NONTARGETING") %>% 
    select(cell, barcode, gene) %>% unique() %>%
    group_by(barcode, gene) %>% tally() %>% filter(n >= 50)

mat = cell_gene %>% group_by(barcode, celltype) %>% 
    tally() %>% dcast(barcode~celltype, fill = 0)
rownames(mat) = as.vector(mat$barcode)
mat = mat[,-1]
mat = as.matrix(mat)

res_NTC = NULL
for(i in 1:nrow(guide_NTC)){
    guide_i = as.vector(guide_NTC$barcode)[i]
    contingency_table = matrix(c(mat[guide_i, ], mat_gene["NONTARGETING", ]), nrow = 2, byrow = TRUE)
    chi_sq_test = chisq.test(contingency_table)
    res_NTC = c(res_NTC, chi_sq_test$p.value)
}
res_NTC = data.frame(barcode = as.vector(guide_NTC$barcode),
                     gene = as.vector(guide_NTC$gene),
                     pval = res_NTC)

res_sig = res[res$pval < quantile(res_NTC$pval, 0.05),]
res_sig = res_sig[order(res_sig$pval),]

saveRDS(list(res, res_NTC, res_sig), paste0(work_path, "/analysis/all_TF_screen_GEx/chi_squared_test_guides_heterogenity.rds"))


### making QQ-plot
res_combine = data.frame(pval = c(res$pval, res_NTC$pval),
                         log10pval = c(-log10(res$pval), -log10(res_NTC$pval)),
                         group = c(rep("TF", nrow(res)), rep("NTC", nrow(res_NTC))))
res_combine$log10pval[is.infinite((res_combine$log10pval))] = 320

qlog10 = function(dat){
    # assuming "nn_de" is a dataframe where each row contains one of your p-values
    n <- length(dat)
    
    # Generate quantiles from a normal distribution
    quantiles <- qnorm(seq(0, 1, length.out = n + 1)[-1])
    
    # Rescale quantiles to fit between 0.00001 and 1
    perfect_values <- pmin(pmax(pnorm(quantiles), 1e-5), 1)
    
    # Apply -log10() function to each value
    log_dis <- -log10(perfect_values)
    
    expected_p <- log_dis[order(-log_dis)]
    
    return(expected_p)
}

dat_1 = subset(res_combine, group == "TF")
dat_1 = dat_1[order(dat_1$pval),]
dat_1$expected_p = qlog10(dat_1$pval)

dat_2 = subset(res_combine, group == "NTC")
dat_2 = dat_2[order(dat_2$pval),]
dat_2$expected_p = qlog10(dat_2$pval)

p = ggplot() + 
    geom_point(data = dat_1, aes(x = expected_p, y = log10pval), color = "blue") + 
    geom_point(data = dat_2, aes(x = expected_p, y = log10pval), color = "red") + 
    geom_hline(yintercept = -log10(quantile(res_NTC$pval, 0.05)), linetype="dotted") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/chi_squared_test_QQplot_guides_heterogenity.pdf"))
print(p)
dev.off()


######################################################################
### Step-9: plot significant examples with multiple significant guides

res_sub = res_sig[res_sig$gene %in% names(table(res_sig$gene)[table(res_sig$gene) > 1]),]
res_sub_2 = head(res_NTC[order(res_NTC$pval),],6)
res_sub_3 = head(res_sig[order(res_sig$pval),],10)
res_sub = rbind(res_sub, res_sub_2, res_sub_3)

for(cnt in 1:nrow(res_sub)){
    print(cnt)
    
    target_i = res_sub$gene[cnt]
    guide_i = res_sub$barcode[cnt]
    
    df = pd %>% select(cell, UMAP_1, UMAP_2) %>% 
        mutate(target = if_else(cell %in% as.vector(cell_gene$cell[cell_gene$barcode == guide_i]), guide_i, "none"))
    
    try(ggplot() +
            geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
            geom_point(data = subset(df, target == guide_i), aes(x = UMAP_1, y = UMAP_2), size=1, color = "red") +
            theme_void() +
            theme(legend.position="none") +
           # labs(title = target_i) +
            theme(plot.title = element_text(hjust = 0.5)) +
            ggsave(paste0(work_path, "/analysis/all_TF_screen_GEx/sgRNA_plot_2/", target_i, "_", guide_i, ".png"),
                   dpi = 300,               
                   height  = 5, 
                   width = 5), silent = TRUE)
}


df = pd %>% select(cell, UMAP_1, UMAP_2) %>% 
    mutate(target = if_else(cell %in% as.vector(cell_gene$cell[cell_gene$gene == "NONTARGETING"]), "NONTARGETING", "none"))

try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(df, target == "NONTARGETING")[sample(1:sum(df$target == "NONTARGETING"), 200),], aes(x = UMAP_1, y = UMAP_2), size=1, color = "red") +
        theme_void() +
        theme(legend.position="none") +
        # labs(title = target_i) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/analysis/all_TF_screen_GEx/sgRNA_plot_2/", "NTC_downsample_200", ".png"),
               dpi = 300,               
               height  = 5, 
               width = 5), silent = TRUE)

