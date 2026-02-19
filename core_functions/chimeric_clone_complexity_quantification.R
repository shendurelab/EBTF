

#######################################################
### Calculate logFC between sgRNA from plasmid and mEBs
### contact: cxqiu@uw.edu

### Loading necessary function and gene meta data

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

#############################
### Step-1: read plasmid file

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tf

plasmid_file = "scaled_plasmid_cut_count.txt"
plasmid_data = read.table(plasmid_file, as.is=T)
plasmid_data$barcode = substr(plasmid_data$gRNA_sequence, 2, nchar(plasmid_data$gRNA_sequence))
plasmid_data = unique(plasmid_data)


white_list = read.table("guide_crispr_cut_scaled_combined_metadata.txt", as.is=T)
plasmid_data_1 = plasmid_data[plasmid_data$barcode %in% white_list$V2,]
plasmid_data_2 = plasmid_data[plasmid_data$gRNA_sequence %in% white_list$V2,]
plasmid_data_2$barcode = as.vector(plasmid_data_2$gRNA_sequence)
plasmid_data = rbind(plasmid_data_1, plasmid_data_2)


########################
### Step-2: read EB file

run_id_list = paste0("EB_", c(1:17,19,21:32))

dat = NULL
for(cnt in 1:length(run_id_list)){
    print(cnt); run_id = run_id_list[cnt]
    dat_i = read.table(paste0(work_path, "/processing/all_TF_screen_sgRNA/", run_id, "_gRNA/cell_gene.txt"),header=T,as.is=T)
    dat = rbind(dat, dat_i)
}

#####################
### Step-3: comparing

plasmid_data$barcode_frac = plasmid_data$guide_count_c_df.plasmid_cut/sum(plasmid_data$guide_count_c_df.plasmid_cut)

df = plasmid_data %>% dplyr::select(barcode, rnames, barcode_frac) %>% rename(plasmid_frac = barcode_frac) %>%
    left_join(dat %>% group_by(barcode) %>% tally(), by = "barcode")
df$n[is.na(df$n)] = 0
df$EB_frac = df$n/sum(df$n)
df$log2_fc = log2(df$EB_frac/df$plasmid_frac)
df = df[!is.na(df$log2_fc),]

print(min(df$log2_fc[!is.infinite((df$log2_fc))]))
df$log2_fc[is.infinite((df$log2_fc))] = -10

p = ggplot(df, aes(log2_fc)) + geom_histogram(bins = 100) +
#    scale_y_continuous(trans='log2') +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf("~/share/log2_fc_EB_plasmid.pdf", 7, 5)
print(p)
dev.off()


#####################################
### Step-4: comparing on target level

plasmid_data_x = plasmid_data %>% group_by(rnames) %>% summarize(count_sum = sum(guide_count_c_df.plasmid_cut))
plasmid_data_x$plasmid_frac = plasmid_data_x$count_sum/sum(plasmid_data_x$count_sum)
plasmid_data_x = plasmid_data_x %>% rename(gene = rnames)

df = plasmid_data_x %>% dplyr::select(gene, plasmid_frac) %>%
    left_join(dat %>% group_by(gene) %>% tally(), by = "gene")
df$n[is.na(df$n)] = 0
df$EB_frac = df$n/sum(df$n)
df$log2_fc = log2(df$EB_frac/df$plasmid_frac)
df = df[!is.na(df$log2_fc),]

print(min(df$log2_fc[!is.infinite((df$log2_fc))]))
df$log2_fc[is.infinite((df$log2_fc))] = -10

white_list = read.table(paste0(work_path, "/processing/all_TF_screen_sgRNA/guide_crispr_cut_scaled_combined_metadata.txt"))

p = ggplot(df, aes(log2_fc)) + geom_histogram(bins = 100) +
    #    scale_y_continuous(trans='log2') +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf("~/share/log2_fc_EB_plasmid_target.pdf", 7, 5)
print(p)
dev.off()

df$plasmid_frac = df$plasmid_frac * 100
df$EB_frac = df$EB_frac * 100
write.table(df[,c("gene", "plasmid_frac", "EB_frac", "log2_fc")], "~/share/log2_fc_EB_plasmid_target.txt", row.names=F, quote=F, sep="\t")


##########################################################
### Step-5: how do the targets in pilot experiment change?

gene_list = c("Foxa1", "Foxa2", "Gata4", "Gata6", "Gbx2", "Greb1l", "Hand1", "Hand2", 
              "Hhex", "Mixl1", "Neurod1", "Otx2", "Pax6", "Prdm1", "Runx1", "Smad1", 
              "Smad2", "Smad3", "Snai2", "Sox1", "Sox17", "T", "Tal1")

df_sub = df %>% filter(gene %in% gene_list) %>% arrange(log2_fc)
summary(2^(df_sub$log2_fc))


###############################################################
### Step-6: can we compare log2(ESC/plasmid) between cut and i?

dat_cut = read.table("guides_table_c.txt", header=T, sep="\t")
dat_cut$target = rownames(dat_cut)
colnames(dat_cut) = c("plasmid", "ES", "target")
dat_cut = dat_cut[,c("target", "plasmid", "ES")]
rownames(dat_cut) = NULL
dat_cut$plasmid_frac = dat_cut$plasmid/sum(dat_cut$plasmid)
dat_cut$ES_frac = dat_cut$ES/sum(dat_cut$ES)
dat_cut$log2_fc = log2(dat_cut$ES_frac/dat_cut$plasmid_frac)
dat_cut$log2_fc[is.infinite(dat_cut$log2_fc)] = -11

x = dat_cut$plasmid[dat_cut$target == "NONTARGETING"]/sum(dat_cut$plasmid)
print(x/(1-x))
x = dat_cut$ES[dat_cut$target == "NONTARGETING"]/sum(dat_cut$ES)
print(x/(1-x))

dat_i = read.table("guide_count_w_guide_names_raw_scaled.txt")
dat_i = dat_i[,c("plasmid_i", "ES_tr_i_scaled", "guide_names")]
dat_i$target = unlist(lapply(as.vector(dat_i$guide_names), function(x) strsplit(x,"[.]")[[1]][1])) 
dat_i = dat_i %>% group_by(target) %>% summarize(plasmid_i_sum = sum(plasmid_i), ES_tr_i_scaled_sum = sum(ES_tr_i_scaled)) %>% as.data.frame()
colnames(dat_i) = c("target", "plasmid", "ES")
dat_i$plasmid_frac = dat_i$plasmid/sum(dat_i$plasmid)
dat_i$ES_frac = dat_i$ES/sum(dat_i$ES)
dat_i$log2_fc = log2(dat_i$ES_frac/dat_i$plasmid_frac)
dat_i$log2_fc[is.infinite(dat_i$log2_fc)] = -11

df = dat_cut %>% select(target, cut_log2_fc = log2_fc) %>%
    full_join(dat_i %>% select(target, i_log2_fc = log2_fc), by = "target")

p = ggplot() + geom_point(data = df, aes(x = cut_log2_fc, y = i_log2_fc)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(x="Log2 fold-change (ESC/plasmid) in CRISPR-cut", y="Log2 fold-change (ESC/plasmid) in CRISPRi", title="") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave("~/share/All_screen_cut_i.pdf", p, width = 5, height = 5)

print(sum(df$cut_log2_fc < 0 & !is.na(df$cut_log2_fc)) / sum(!is.na(df$cut_log2_fc)))
print(sum(df$i_log2_fc < 0 & !is.na(df$i_log2_fc)) / sum(!is.na(df$i_log2_fc)))

df_o = df
df_o$cut_log2_fc = round(df$cut_log2_fc, 2)
df_o$i_log2_fc = round(df$i_log2_fc, 2)
write.table(df_o, "~/share/All_screen_cut_i.txt", row.names=F, sep="\t", quote=F)



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







