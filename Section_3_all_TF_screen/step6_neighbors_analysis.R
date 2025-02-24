
#############################################################
### Performing lochness analysis on the big screen experiment
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
    return(mutant_neighbors / kadj)
}


#########################################################
### Step-1: calculating neighbors scores for each targets

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

obj = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned.rds"))
pd = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned_pd.rds"))
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
saveRDS(nn_matrix, paste0(work_path, "/analysis/all_TF_screen_GEx/neighbors_result/nn_matrix_k20.rds"))

cell_gene = readRDS(paste0(work_path, "/processing/all_TF_screen_sgRNA/cell_gene.rds"))
cell_gene$cell = paste0("all_TF_screen_GEx_", as.vector(cell_gene$cell))

cell_gene_include = cell_gene %>% group_by(gene) %>% tally() %>% filter(n >= 50)

gene_list = unique(cell_gene_include$gene)
gene_list = gene_list[gene_list != "NONTARGETING"]
### n = 820 targets

res = NULL
for(i in 1:length(gene_list)){
    gene_i = gene_list[i]
    print(paste0(i, "/", length(gene_list), " : ", gene_i))
    
    cell_include_i = as.vector(cell_gene$cell[cell_gene$gene == gene_i])
    cell_include_i = c(1:nrow(nn_matrix))[rownames(nn_matrix) %in% cell_include_i]
    
    ### downsampling to 200 cells across targets
    if(length(cell_include_i) > 200){
        cell_include_i = sample(cell_include_i, 200)
    }
   
    nn_matrix_i = nn_matrix[cell_include_i,]
    
    res_i = apply(nn_matrix_i, 1, function(x) mean(x %in% cell_include_i))
    res = c(res, mean(res_i))
}

res = data.frame(gene = gene_list,
                 score = res)

saveRDS(res, paste0(work_path, "/analysis/all_TF_screen_GEx/neighbors_result/res_neighbors_TF.rds"))


### comparing the result with chi-squared test

y = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/chi_squared_test.rds"))
y = y[[1]]
y$log10pval = -log10(y$pval)
y$log10pval[is.infinite(y$log10pval)] = 320

cor.test(res$score, y$log10pval, method = "spearman")
### rho = 0.588734; p-value = 1.204152e-77

df = data.frame(gene = res$gene, score = res$score, log10pval = y$log10pval)
p = ggplot() + 
    geom_point(data = df, aes(x = score, y = log10pval)) + 
    geom_point(data = subset(df, log10pval > (-log10(1.26e-32))),
               aes(x = score, y = log10pval), color = "blue") + 
    geom_text(data=subset(df, score > 0.05 & log10pval > (-log10(1.26e-32))),
              aes(score, log10pval, label=gene), hjust = 0.75, color = "red", size = 2) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/comparing_chisquared_neighbors.pdf"))
print(p)
dev.off()






################################################################################
### Step-2: repeating the neighbors score, by randomly selecting three guides of NTC, and repeat it 5000 times


obj = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned.rds"))
pd = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned_pd.rds"))

cell_gene = readRDS(paste0(work_path, "/processing/all_TF_screen_sgRNA/cell_gene.rds"))
cell_gene$cell = paste0("all_TF_screen_GEx_", as.vector(cell_gene$cell))

NTC_guide_num = cell_gene %>% filter(gene == "NONTARGETING") %>% 
    group_by(barcode) %>% tally() %>% filter(n >= 20)
### n = 540 guides

### how many neighbors are considering 
kadj = 20

res = NULL

iter_time = 1

for(cnt in 1:iter_time){
    print(paste0(cnt, "/", iter_time))
    
    each_group = 3
    group_num = floor(nrow(NTC_guide_num)/each_group) 
    set.seed(1234)
    NTC_guide_num$group = sample(rep(1:group_num, each = each_group))[1:nrow(NTC_guide_num)]
    
    for(group_i in 1:group_num){
        print(paste0(group_i, "/", group_num))
        guides_include = as.vector(NTC_guide_num$barcode[NTC_guide_num$group == group_i])
        cell_include = as.vector(cell_gene$cell[cell_gene$barcode %in% guides_include])
        if(length(cell_include) > 100){
            cell_include = cell_include[sample(1:length(cell_include), 100)]
        }
        obj$genotype = ifelse(rownames(pd) %in% cell_include, "KO", "WT")
        res = cbind(res, calculate_neighbors(obj, kadj))
    }
}

saveRDS(res, paste0(work_path, "/analysis/all_TF_screen_GEx/neighbors_result/neighbors_NTC.rds"))





