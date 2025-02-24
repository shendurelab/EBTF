
###############################
### Performing chi-squared test
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tf
### validated_125TF.cell_gene.rds

pd = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned_pd.rds"))
pd$cell = gsub("EB21_CRISPRcut_125TF_GEx_", "", as.vector(pd$cell_id))
cell_gene = readRDS("validated_125TF.cell_gene.rds")

####################################################################################
### Step-1: regroup cells, to ensure more cells in each cell cluster (not necessary)
celltype_regroup = as.vector(pd$celltype)
celltype_regroup[pd$celltype %in% c("Neuroectoderm",
                                    "Early neurons",
                                    "Floor plate",
                                    "Eye field")] = "Neuroectoderm"
pd$celltype_regroup = as.vector(celltype_regroup) 

pd$replicate = rep("rep1", nrow(pd))
pd$replicate[pd$experiment_id %in% paste0("EB21_CRISPRcut_125TF_GEx_EB_", c(17:32))] = "rep2"

cell_gene = cell_gene %>% left_join(pd[,c("cell","celltype","celltype_regroup","replicate")], by = "cell")


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
print(nrow(target_cell_num)) ### n = 118

### filtering out NTCs with less than 20 cells

NTC_guide_num = NTC_guide_num %>% filter(n >= 20)
print(nrow(NTC_guide_num)) ### n = 80




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
res = data.frame(gene = as.vector(target_cell_num$gene),
                 pval = res)


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
    group_num = floor(nrow(NTC_guide_num)/each_group) + 1
    set.seed(1234)
    NTC_guide_num$group = sample(rep(1:group_num, each = each_group))[1:nrow(NTC_guide_num)]
    
    for(i in 1:group_num){
        barcode_include = as.vector(NTC_guide_num$barcode[NTC_guide_num$group == i])
        
        mat_sub = cell_gene %>% filter(barcode %in% barcode_include) %>% 
#            slice_sample(n = 100) %>%
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
saveRDS(list(res, res_NTC, res_sig), paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/chi_squared_test.rds"))


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
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/chi_squared_test_QQplot_x.pdf"))
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
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/chi_squared_test_QQplot.pdf"))
print(p)
dev.off()



###########################################
### Step-6: plot gRNA for individual target


pd = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned_pd.rds"))
cell_gene = readRDS("validated_125TF.cell_gene.rds")
pd$cell = gsub("EB21_CRISPRcut_125TF_GEx_", "", as.vector(pd$cell_id))

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
            ggsave(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/sgRNA_plot/", target_i, ".png"),
                   dpi = 300,               
                   height  = 5, 
                   width = 5), silent = TRUE)
}




###################################################################
### Step-7: specifically, plot Carm1 gRNAs distribution in the UMAP

guide_i = "Carm1"

df = pd %>% select(cell, UMAP_1, UMAP_2, replicate) %>% 
    mutate(target = if_else(cell %in% as.vector(cell_gene$cell[cell_gene$gene == guide_i]), guide_i, "none"))

try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(df, target == guide_i & replicate == "rep1"), aes(x = UMAP_1, y = UMAP_2), size=0.3, color = "red") +
        theme_void() +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/Carm1_gRNA_umap_rep1.png"),
               dpi = 300,               
               height  = 3.7, 
               width = 4), silent = TRUE)


try(ggplot() +
        geom_point(data = df, aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
        geom_point(data = subset(df, target == guide_i & replicate == "rep2"), aes(x = UMAP_1, y = UMAP_2), size=0.3, color = "blue") +
        theme_void() +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/Carm1_gRNA_umap_rep2.png"),
               dpi = 300,               
               height  = 3.7, 
               width = 4), silent = TRUE)





##########################################
### Step-8: comparing to all_TF_screen_GEx

dat_all = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/chi_squared_test.rds"))
dat_125TF = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/chi_squared_test.rds"))

res_all = dat_all[[1]]
res_125TF = dat_125TF[[1]]

dat_combine = res_all %>% rename(pval_all = pval) %>%
    left_join(res_125TF %>% rename(pval_125TF = pval), by = "gene") %>%
    filter(!is.na(pval_125TF))

dat_combine$log10_pval_all = -log10(dat_combine$pval_all)
dat_combine$log10_pval_125TF = -log10(dat_combine$pval_125TF)

print(max(dat_combine$log10_pval_all[!is.infinite((dat_combine$log10_pval_all))]))
print(max(dat_combine$log10_pval_125TF[!is.infinite((dat_combine$log10_pval_125TF))]))

dat_combine$log10_pval_all[is.infinite((dat_combine$log10_pval_all))] = 320
dat_combine$log10_pval_125TF[is.infinite((dat_combine$log10_pval_125TF))] = 320

cor.test(dat_combine$log10_pval_all, dat_combine$log10_pval_125TF, method = "spearman")

p = ggplot() + 
    geom_point(data = dat_combine, aes(x = log10_pval_all, y = log10_pval_125TF)) + 
    #geom_point(data = subset(dat_combine, pval_all < quantile(dat_all[[2]], 0.05)), aes(x = log10_pval_all, y = log10_pval_125TF), color = "red") + 
    #geom_point(data = subset(dat_combine, pval_125TF < quantile(dat_125TF[[2]], 0.05)), aes(x = log10_pval_all, y = log10_pval_125TF), color = "blue") + 
    geom_text(data = subset(dat_combine, log10_pval_all > quantile(dat_combine$log10_pval_all, 0.9)),
              aes(log10_pval_all, log10_pval_125TF, label=gene), hjust = 0.75, color = "red", size = 2) +
    geom_text(data = subset(dat_combine, log10_pval_125TF > quantile(dat_combine$log10_pval_125TF, 0.9)),
              aes(log10_pval_all, log10_pval_125TF, label=gene), hjust = 0.75, color = "blue", size = 2) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/comparing_two_datasets.pdf"))
print(p)
dev.off()


dat_combine_two_datasets = dat_combine



####################################
### Step-9: comparing two replicates

target_cell_num_rep1 = cell_gene %>% filter(gene != "NONTARGETING", replicate == "rep1") %>% 
    select(cell, gene) %>% unique() %>%
    group_by(gene) %>% tally() %>% filter(n >= 50)

target_cell_num_rep2 = cell_gene %>% filter(gene != "NONTARGETING", replicate == "rep2") %>% 
    select(cell, gene) %>% unique() %>%
    group_by(gene) %>% tally() %>% filter(n >= 50)

mat = cell_gene %>% filter(replicate == "rep1") %>% group_by(gene, celltype) %>% 
    tally() %>% dcast(gene~celltype, fill = 0)
rownames(mat) = as.vector(mat$gene)
mat = mat[,-1]
mat = as.matrix(mat)

mat_rep1 = mat

res = NULL
for(i in 1:nrow(target_cell_num_rep1)){
    target_i = as.vector(target_cell_num_rep1$gene)[i]
    contingency_table = matrix(c(mat[target_i, ], mat["NONTARGETING", ]), nrow = 2, byrow = TRUE)
    chi_sq_test = chisq.test(contingency_table)
    res = c(res, chi_sq_test$p.value)
}
res = data.frame(gene = as.vector(target_cell_num_rep1$gene),
                 pval = res)
res_rep1 = res

mat = cell_gene %>% filter(replicate == "rep2") %>% group_by(gene, celltype) %>% 
    tally() %>% dcast(gene~celltype, fill = 0)
rownames(mat) = as.vector(mat$gene)
mat = mat[,-1]
mat = as.matrix(mat)

mat_rep2 = mat

res = NULL
for(i in 1:nrow(target_cell_num_rep2)){
    target_i = as.vector(target_cell_num_rep2$gene)[i]
    contingency_table = matrix(c(mat[target_i, ], mat["NONTARGETING", ]), nrow = 2, byrow = TRUE)
    chi_sq_test = chisq.test(contingency_table)
    res = c(res, chi_sq_test$p.value)
}
res = data.frame(gene = as.vector(target_cell_num_rep2$gene),
                 pval = res)
res_rep2 = res


dat_combine = res_rep1 %>% rename(pval_rep1 = pval) %>%
    left_join(res_rep2 %>% rename(pval_rep2 = pval), by = "gene") %>%
    filter(!is.na(pval_rep2))

dat_combine$log10_pval_rep1 = -log10(dat_combine$pval_rep1)
dat_combine$log10_pval_rep2 = -log10(dat_combine$pval_rep2)

print(max(dat_combine$log10_pval_rep1[!is.infinite((dat_combine$log10_pval_rep1))]))
print(max(dat_combine$log10_pval_rep2[!is.infinite((dat_combine$log10_pval_rep2))]))

dat_combine$log10_pval_rep1[is.infinite((dat_combine$log10_pval_rep1))] = 165
dat_combine$log10_pval_rep2[is.infinite((dat_combine$log10_pval_rep2))] = 90

cor.test(dat_combine$log10_pval_rep1, dat_combine$log10_pval_rep2, method = "spearman")

p = ggplot() + 
    geom_point(data = dat_combine, aes(x = log10_pval_rep1, y = log10_pval_rep2)) + 
    geom_text(data = subset(dat_combine, log10_pval_rep1 > quantile(dat_combine$log10_pval_rep1, 0.9)),
              aes(log10_pval_rep1, log10_pval_rep2, label=gene), hjust = 0.75, color = "red", size = 2) +
    geom_text(data = subset(dat_combine, log10_pval_rep2 > quantile(dat_combine$log10_pval_rep2, 0.9)),
              aes(log10_pval_rep1, log10_pval_rep2, label=gene), hjust = 0.75, color = "blue", size = 2) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/comparing_two_replicates.pdf"))
print(p)
dev.off()


dat_combine_two_replicates = dat_combine



#####################################################################################################
### Step-10: calculating odds ratio for the top TF targets for individual cell, across three datasets

### validated dataset - rep1
mat_rep1

mat_rep1_neuroectoderm = apply(mat_rep1[,c("Neuroectoderm",
                                           "Early neurons",
                                           "Floor plate",
                                           "Eye field")], 1, sum)
mat_rep1 = mat_rep1[,!colnames(mat_rep1) %in% c("Neuroectoderm",
                                                "Early neurons",
                                                "Floor plate",
                                                "Eye field")]
mat_rep1 = cbind(mat_rep1, mat_rep1_neuroectoderm)
colnames(mat_rep1) = c(colnames(mat_rep1)[1:15], "Neuroectoderm")

### validated dataset - rep2
mat_rep2

mat_rep2_neuroectoderm = apply(mat_rep2[,c("Neuroectoderm",
                                           "Early neurons",
                                           "Floor plate",
                                           "Eye field")], 1, sum)
mat_rep2 = mat_rep2[,!colnames(mat_rep2) %in% c("Neuroectoderm",
                                                "Early neurons",
                                                "Floor plate",
                                                "Eye field")]
mat_rep2 = cbind(mat_rep2, mat_rep2_neuroectoderm)
colnames(mat_rep2) = c(colnames(mat_rep2)[1:15], "Neuroectoderm")


### big screen dataset
pd = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned_pd.rds"))
pd$cell = gsub("all_TF_screen_GEx_", "", as.vector(pd$cell_id))
cell_gene = readRDS(paste0(work_path, "/processing/all_TF_screen_sgRNA/cell_gene.rds"))
cell_gene = cell_gene %>% left_join(pd[,c("cell","celltype")], by = "cell")

mat = cell_gene %>% group_by(gene, celltype) %>% 
    tally() %>% dcast(gene~celltype, fill = 0)
rownames(mat) = as.vector(mat$gene)
mat = mat[,-1]
mat = as.matrix(mat)

mat_dat1 = mat


### which TFs are presented in the heatmaps (top 3 genes in each dataset)
gene_include = intersect(dat_combine_two_replicates$gene, dat_combine_two_datasets$gene)
TF_1 = dat_combine_two_datasets %>% filter(gene %in% gene_include) %>%
    slice_max(order_by = log10_pval_all, n = 5)
TF_2 = dat_combine_two_replicates %>% filter(gene %in% gene_include) %>%
    slice_max(order_by = log10_pval_rep1, n = 5)
TF_3 = dat_combine_two_replicates %>% filter(gene %in% gene_include) %>%
    slice_max(order_by = log10_pval_rep2, n = 5)

TF_include = c("Elf2", "Rfx2", "Zfp958", "Carm1", "Batf2", "Dlx4", "Hoxa4", "Elk3", "Dmrt2")

### calculating odds ratio in each dataset for those TF included
### big screen dataset
mat = mat_dat1
mat = mat[c(TF_include, "NONTARGETING"),]
res_or = NULL
res_pval = NULL
for(i in TF_include){
    or_i = NULL
    pval_i = NULL
    mat_1 = mat[i,]; mat_1_sum = sum(mat_1)
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
rownames(res_or) = rownames(res_pval) = paste0("dat1_", TF_include)
colnames(res_or) = colnames(res_pval) = colnames(mat)
res_or_dat1 = res_or


### validated dataset - rep1
mat = mat_rep1
mat = mat[c(TF_include, "NONTARGETING"),]
res_or = NULL
res_pval = NULL
for(i in TF_include){
    or_i = NULL
    pval_i = NULL
    mat_1 = mat[i,]; mat_1_sum = sum(mat_1)
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
rownames(res_or) = rownames(res_pval) = paste0("rep1_", TF_include)
colnames(res_or) = colnames(res_pval) = colnames(mat)
res_or_rep1 = res_or

### validated dataset - rep2
mat = mat_rep2
mat = mat[c(TF_include, "NONTARGETING"),]
res_or = NULL
res_pval = NULL
for(i in TF_include){
    or_i = NULL
    pval_i = NULL
    mat_1 = mat[i,]; mat_1_sum = sum(mat_1)
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
rownames(res_or) = rownames(res_pval) = paste0("rep2_", TF_include)
colnames(res_or) = colnames(res_pval) = colnames(mat)
res_or_rep2 = res_or


### combining three odds ratio matrix and making the heatmap
col_names = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% colnames(res_or_dat1)]
aggr_matrix = rbind(res_or_dat1[,col_names],
                    res_or_rep1[,col_names],
                    res_or_rep2[,col_names])
row_names = paste0(rep(c("dat1_", "rep1_", "rep2_"), times = length(TF_include)), rep(TF_include, each = 3))
aggr_matrix = aggr_matrix[row_names,]


library("gplots")
library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/Heatmap_three_datasets.pdf"), 8, 5)
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





