
##################################################################
### Performing lochness analysis on the validated 125TF experiment
### Chengxiang Qiu
### Feb-20, 2025

calculate_lochness = function(obj, kadj){
    
    safe_column_multiply = function(bmat, scale_vector) {
        bmat@x <- bmat@x * rep.int(scale_vector, diff(bmat@p))
        return(bmat)
    }
    # kadj = round(0.5 * sqrt(ncol(obj)))
    
    # MT scores
    mt_counts = as.data.frame(table(obj$genotype))
    mutant_mask = ifelse(obj@meta.data$genotype=="WT", 0, 1)
    mutant_neighbors = Matrix::rowSums(safe_column_multiply(obj@graphs$RNA_nn, mutant_mask))
    pca_mutant_score = (mutant_neighbors / kadj) / (1 - mt_counts$Freq[mt_counts$Var1=='WT'] / length(mutant_neighbors)) - 1
    
    return(pca_mutant_score)
}


########################################################
### Step-1: calculating lochness scores for each targets

obj = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned.rds"))
pd = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned_pd.rds"))
pca_coor = Embeddings(obj, reduction = "pca")

### how many neighbors are considering 
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
saveRDS(nn_matrix, paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/neighbors_result/nn_matrix_k20.rds"))



cell_gene = readRDS(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/cell_gene.rds"))
cell_gene$cell = paste0("EB21_CRISPRcut_125TF_GEx_", as.vector(cell_gene$cell))

cell_gene_include = cell_gene %>% group_by(gene) %>% tally() %>% filter(n >= 50)

gene_list = unique(cell_gene_include$gene)
gene_list = gene_list[gene_list != "NONTARGETING"]
### n = 118 targets

res = NULL
for(i in 1:length(gene_list)){
    gene_i = gene_list[i]
    print(paste0(i, "/", length(gene_list), " : ", gene_i))
    
    cell_include_i = as.vector(cell_gene$cell[cell_gene$gene == gene_i])
    cell_include_i = c(1:nrow(nn_matrix))[rownames(nn_matrix) %in% cell_include_i]
    
    ### downsampling to 1000 cells across targets
    if(length(cell_include_i) > 1000){
        cell_include_i = sample(cell_include_i, 1000)
    }
    
    nn_matrix_i = nn_matrix[cell_include_i,]
    
    res_i = apply(nn_matrix_i, 1, function(x) mean(x %in% cell_include_i))
    res = c(res, mean(res_i))
}

res = data.frame(gene = gene_list,
                 score = res)

saveRDS(res, paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/neighbors_result/res_neighbors_TF.rds"))


######################################################
### Step-2: comparing the result with chi-squared test

y = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/chi_squared_test.rds"))
y = y[[1]]
y$log10pval = -log10(y$pval)
y$log10pval[is.infinite(y$log10pval)] = 320

cor.test(res$score, y$log10pval, method = "spearman")
### rho = 0.6826421; p-value = 1.689838e-17

df = data.frame(gene = res$gene, score = res$score, log10pval = y$log10pval)
p = df %>%
    ggplot(aes(x = score, y = log10pval)) + 
    geom_point() + 
    geom_text(data=subset(df, score > 0.05 & log10pval > quantile(y$log10pval, 0.95)),
              aes(score, log10pval, label=gene), hjust = 0.75, color = "red", size = 2) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/comparing_chisquared_neighbors.pdf"))
print(p)
dev.off()





##########################################
### Step-3: comparing to all_TF_screen_GEx

res_all = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/neighbors_result/res_neighbors_TF.rds"))
res_125TF = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/neighbors_result/res_neighbors_TF.rds"))


dat_combine = res_all %>% rename(score_all = score) %>%
    left_join(res_125TF %>% rename(score_125TF = score), by = "gene") %>%
    filter(!is.na(score_125TF))

cor.test(dat_combine$score_all, dat_combine$score_125TF, method = "spearman")

p = ggplot() + 
    geom_point(data = dat_combine, aes(x = score_all, y = score_125TF)) + 
    geom_text(data = subset(dat_combine, score_all > quantile(dat_combine$score_all, 0.9)),
              aes(score_all, score_125TF, label=gene), hjust = 0.75, color = "red", size = 2) +
    geom_text(data = subset(dat_combine, score_125TF > quantile(dat_combine$score_125TF, 0.9)),
              aes(score_all, score_125TF, label=gene), hjust = 0.75, color = "blue", size = 2) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/comparing_two_datasets_neighbors_analysis.pdf"))
print(p)
dev.off()







####################################
### Step-4: comparing two replicates



pd = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned_pd.rds"))
pd$cell = gsub("EB21_CRISPRcut_125TF_GEx_", "", as.vector(pd$cell_id))
cell_gene = readRDS(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/cell_gene.rds"))

nn_matrix = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/neighbors_result/nn_matrix_k20.rds"))
rownames(nn_matrix) = as.vector(pd$cell)

pd$replicate = rep("rep1", nrow(pd))
pd$replicate[pd$experiment_id %in% paste0("EB21_CRISPRcut_125TF_GEx_EB_", c(17:32))] = "rep2"
cell_gene = cell_gene %>% left_join(pd[,c("cell","celltype","replicate")], by = "cell")

target_cell_num_rep1 = cell_gene %>% filter(gene != "NONTARGETING", replicate == "rep1") %>% 
    select(cell, gene) %>% unique() %>%
    group_by(gene) %>% tally() %>% filter(n >= 50)

target_cell_num_rep2 = cell_gene %>% filter(gene != "NONTARGETING", replicate == "rep2") %>% 
    select(cell, gene) %>% unique() %>%
    group_by(gene) %>% tally() %>% filter(n >= 50)

gene_list = unique(target_cell_num_rep1$gene)
res = NULL
for(i in 1:length(gene_list)){
    gene_i = gene_list[i]
    print(paste0(i, "/", length(gene_list), " : ", gene_i))
    
    cell_include_i = as.vector(cell_gene$cell[cell_gene$gene == gene_i & cell_gene$replicate == "rep1"])
    cell_include_i = c(1:nrow(nn_matrix))[rownames(nn_matrix) %in% cell_include_i]
    
    ### downsampling to 500 cells across targets
    if(length(cell_include_i) > 500){
        cell_include_i = sample(cell_include_i, 500)
    }
    
    
    nn_matrix_i = nn_matrix[cell_include_i,]
    
    res_i = apply(nn_matrix_i, 1, function(x) mean(x %in% cell_include_i))
    res = c(res, mean(res_i))
}
res_rep1 = data.frame(gene = gene_list,
                      score = res)

gene_list = unique(target_cell_num_rep2$gene)
res = NULL
for(i in 1:length(gene_list)){
    gene_i = gene_list[i]
    print(paste0(i, "/", length(gene_list), " : ", gene_i))
    
    cell_include_i = as.vector(cell_gene$cell[cell_gene$gene == gene_i & cell_gene$replicate == "rep2"])
    cell_include_i = c(1:nrow(nn_matrix))[rownames(nn_matrix) %in% cell_include_i]
    
    ### downsampling to 500 cells across targets
    if(length(cell_include_i) > 500){
        cell_include_i = sample(cell_include_i, 500)
    }
    
    nn_matrix_i = nn_matrix[cell_include_i,]
    
    res_i = apply(nn_matrix_i, 1, function(x) mean(x %in% cell_include_i))
    res = c(res, mean(res_i))
}
res_rep2 = data.frame(gene = gene_list,
                      score = res)


dat_combine = res_rep1 %>% rename(score_rep1 = score) %>%
    left_join(res_rep2 %>% rename(score_rep2 = score), by = "gene") %>%
    filter(!is.na(score_rep2))

cor.test(dat_combine$score_rep1, dat_combine$score_rep2, method = "spearman")

p = ggplot() + 
    geom_point(data = dat_combine, aes(x = score_rep1, y = score_rep2)) + 
    geom_text(data = subset(dat_combine, score_rep1 > quantile(dat_combine$score_rep1, 0.9)),
              aes(score_rep1, score_rep2, label=gene), hjust = 0.75, color = "red", size = 2) +
    geom_text(data = subset(dat_combine, score_rep2 > quantile(dat_combine$score_rep2, 0.9)),
              aes(score_rep1, score_rep2, label=gene), hjust = 0.75, color = "blue", size = 2) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/comparing_two_replicates_neighbors_analysis.pdf"))
print(p)
dev.off()








