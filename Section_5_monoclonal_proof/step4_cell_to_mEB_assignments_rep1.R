
####################################################
### Assiging cells to each mEBs based on the barcode
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tf
### NTC_scaled_sgRNA_BC_final_list_20210510.txt

white_list = read.table("NTC_scaled_sgRNA_BC_final_list_20210510.txt", header=T, as.is=T)
white_list$sgRNA = paste0("G", substr(white_list$sgRNA, 1, nchar(white_list$sgRNA)-1))

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/obj_processed_pd.rds"))

run_list = paste0("E2E_rep", c("1_1", "1_2"))

################################################
### Step-1: extracting and filtering BC and gRNA 

bc_matrix = NULL

for(run_id in run_list){    
    print(run_id)
    pd_filter = subset(pd, experiment_id == paste0("monoclonal_EB_proof_", run_id))
    
    bc_matrix_table = read.table(paste0(data_path, "/data/monoclonal_EB_proof/", run_id, "_cell_bc.txt"), header=T, as.is=T)
    bc_matrix_table$cell_id = paste0("monoclonal_EB_proof_", run_id, "_",  unlist(lapply(as.vector(bc_matrix_table$cBC), function(x) strsplit(x,"[-]")[[1]][1])))
    bc_matrix_table = subset(bc_matrix_table, cell_id %in% as.vector(pd_filter$cell_id))
    
    bc_matrix = rbind(bc_matrix, bc_matrix_table)
}

tmp_rep = bc_matrix %>% group_by(mBC) %>%
    summarize(read_count_sum = sum(n_reads_filtered), umi_count_sum = sum(filtered_corrected_UMIs)) %>%
    mutate(log2_read_count = log2(read_count_sum), log2_umi_count = log2(umi_count_sum))

try(ggplot() +
        geom_point(data = subset(tmp_rep, log2_read_count < 15), aes(x = log2_read_count, y = log2_umi_count), size = 0.5) + 
        geom_point(data = subset(tmp_rep, log2_read_count >= 15), aes(x = log2_read_count, y = log2_umi_count), color = "red", size = 1) +
        geom_vline(xintercept = 15) +
        labs(x="Log2 read count per BC", y="Log2 UMI count per BC", title="") +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/BC_2D_hist_rep1.pdf"),
               height  = 3, 
               width = 5), silent = TRUE)

try(ggplot() +
        geom_point(data = subset(tmp_rep, log2_read_count < 15), aes(x = log2_read_count, y = log2_umi_count), size = 1) + 
        geom_point(data = subset(tmp_rep, log2_read_count >= 15), aes(x = log2_read_count, y = log2_umi_count), color = "red", size = 2) +
        labs(x="Log2 read count per BC", y="Log2 UMI count per BC", title="") +
        theme_void() +
        theme(legend.position="none") +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/BC_2D_hist_rep1.png"),
               dpi = 300,
               height  = 3, 
               width = 5), silent = TRUE)

### After checking the 2D-hist, I decided to use log2_read_count >= 15 as cutoff to filter BCs
### n = 26 BCs are retained
tmp_include = subset(tmp_rep, log2_read_count >= 15)

### filtering BCs and gRNAs
bc_matrix_list = list()
gRNA_matrix_list = list()

for(run_id in run_list){    
    print(run_id)
    pd_filter = subset(pd, experiment_id == paste0("monoclonal_EB_proof_", run_id))
    
    count_matrix_all = Read10X(paste0(data_path, "/data/monoclonal_EB_proof/", run_id, "/outs/raw_feature_bc_matrix"), 
                               gene.column = 1,
                               strip.suffix = T)
    gRNA_matrix = count_matrix_all[["CRISPR Guide Capture"]]
    colnames(gRNA_matrix) = paste0("monoclonal_EB_proof_", run_id, "_", colnames(gRNA_matrix))
    gRNA_matrix = gRNA_matrix[,colnames(gRNA_matrix) %in% as.vector(pd_filter$cell_id)]
    gRNA_matrix = gRNA_matrix[rowMeans(gRNA_matrix > 0) >= 0.01,]
    
    gRNA_matrix_list[[run_id]] = gRNA_matrix
    
    bc_matrix_table = read.table(paste0(data_path, "/data/monoclonal_EB_proof/", run_id, "_cell_bc.txt"), header=T, as.is=T)
    bc_matrix_table$cell_id = paste0("monoclonal_EB_proof_", run_id, "_",  unlist(lapply(as.vector(bc_matrix_table$cBC), function(x) strsplit(x,"[-]")[[1]][1])))
    bc_matrix_table = subset(bc_matrix_table, cell_id %in% as.vector(pd_filter$cell_id) & mBC %in% as.vector(tmp_include$mBC))
    
    bc_matrix = bc_matrix_table[,c("mBC", "cell_id", "filtered_corrected_UMIs")]
    bc_matrix = dcast(bc_matrix, mBC ~ cell_id, fill = 0)
    rownames(bc_matrix) = bc_matrix[,1]
    bc_matrix = bc_matrix[,-1]
    bc_matrix = as(bc_matrix, "sparseMatrix")
    
    bc_matrix_list[[run_id]] = bc_matrix
}

### double-checking if feature names are matched between replicates
### n = 26 BCs, n = 22 gRNAs
print(sum(rownames(bc_matrix_list[[1]]) %in% rownames(bc_matrix_list[[2]])))
print(sum(rownames(gRNA_matrix_list[[1]]) %in% rownames(gRNA_matrix_list[[2]])))

bc_matrix = NULL
gRNA_matrix = NULL

for(run_id in run_list){
    if(is.null(bc_matrix)){
        bc_matrix = bc_matrix_list[[run_id]]
    } else {
        bc_matrix = cbind(bc_matrix, bc_matrix_list[[run_id]][rownames(bc_matrix),])
    }
    
    if(is.null(gRNA_matrix)){
        gRNA_matrix = gRNA_matrix_list[[run_id]]
    } else {
        gRNA_matrix = cbind(gRNA_matrix, gRNA_matrix_list[[run_id]][rownames(gRNA_matrix),])
    }
}

saveRDS(bc_matrix, paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/BC_matrix_rep1.rds"))
saveRDS(gRNA_matrix, paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/gRNA_matrix_rep1.rds"))



#########################################################################
### Step-2: Setting cutoffs to assign cells to clonotype based on the BCs

cell_BC = 100*t(t(bc_matrix)/colSums(bc_matrix))
cell_BC = melt(as.matrix(cell_BC))
names(cell_BC) = c("BC", "cell_id", "pct")
cell_BC_x = melt(as.matrix(bc_matrix))
cell_BC$UMI_count = as.vector(cell_BC_x$value)
cell_BC = cell_BC %>% group_by(cell_id) %>% arrange(desc(pct), .by_group = TRUE)
cell_BC$pct_order = rep(c(1:nrow(bc_matrix)), times = ncol(bc_matrix))

ggplot(data = subset(cell_BC, UMI_count > 0), aes(log2(UMI_count))) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = log2(10)) + theme_classic(base_size = 10) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/hist_umi_per_BC_per_cell_rep1.pdf"),
           height  = 4, 
           width = 4)

ggplot(data = subset(cell_BC, UMI_count > 0), aes(pct)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = 30) + theme_classic(base_size = 10) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/hist_pct_umi_per_BC_per_cell_rep1.pdf"),
           height  = 4, 
           width = 4)


### pct >= 30 and UMI_count >= 10
cell_BC_sub = subset(cell_BC, pct >= 30 & UMI_count >= 10)
print(length(unique(cell_BC_sub$cell_id))/length(unique(cell_BC$cell_id)))
print(sum(table(cell_BC_sub$cell_id) == 1)/length(unique(cell_BC$cell_id)))
### 86.9% cells are retained, and 83.0% cells have unique BC




##############################################################
### Step-3: calculating Peasrson correlation on the UMI matrix

BC_list_num = cell_BC_sub %>% group_by(BC) %>% tally() %>% arrange(desc(n)) %>% 
    mutate(log2_n = log2(n)) %>% as.data.frame()
BC_list_num$BC_id = paste0("BC_",1:nrow(BC_list_num))
BC_list = as.vector(BC_list_num$BC)

BC_list_num$BC = factor(BC_list_num$BC, levels = BC_list)
p = ggplot(data=BC_list_num, aes(x=BC, y=n)) +
    geom_bar(stat="identity") + theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_fill_viridis() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    #scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = "cell_count")) +
    #guides(y = "none") +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/BC_cell_num.pdf"),
           height  = 5, 
           width = 5)


### calculating Peasrson correlation on the UMI matrix
bc_matrix_sub = bc_matrix[,colnames(bc_matrix) %in% cell_BC_sub$cell_id]
bc_matrix_sub_norm = t(bc_matrix_sub)/colSums(bc_matrix_sub)
bc_matrix_sub_scale = scale(bc_matrix_sub_norm)
bc_matrix_sub_scale = bc_matrix_sub_scale[,BC_list]

dat = NULL
for(i in 1:length(BC_list)){
    print(i)
    for(j in 1:length(BC_list)){
        dat = rbind(dat, data.frame(A = i,
                                    B = j,
                                    corr = cor.test(bc_matrix_sub_scale[,i], bc_matrix_sub_scale[,j])$estimate))
    }
}

dat = subset(dat, A >= B)

p = ggplot() + 
    geom_point(data = subset(dat, A == B), aes(A, B), size=6, shape=22, color="grey80", fill="white") + 
    geom_point(data = subset(dat, A != B), aes(A, B, fill = corr), size=6, shape=22, color="grey80") + 
    labs(x="", y="", title="") +
    theme_classic(base_size = 10) +
    #theme(legend.position="none") +
    scale_fill_viridis() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/heatmap_UMI_correlation_BC.pdf"),
           height  = 5, 
           width = 6)

### identifying clones with MOI > 1
dat %>% filter(A != B, corr > 0) %>% arrange(desc(corr)) 

###########################################
### Step-4: identifying distinct clonotypes

cell_BC_sub_uniq_x = cell_BC_sub %>% group_by(cell_id) %>% tally() %>% filter(n == 1)

cell_BC_sub_uniq = cell_BC_sub %>% filter(cell_id %in% as.vector(cell_BC_sub_uniq_x$cell_id)) %>%
    left_join(BC_list_num[,c("BC","BC_id")], by = "BC")

cell_BC_sub_mult = cell_BC_sub %>% filter(!cell_id %in% as.vector(cell_BC_sub_uniq_x$cell_id)) %>%
    left_join(BC_list_num[,c("BC","BC_id")], by = "BC")
cell_BC_sub_mult$BC_id = factor(cell_BC_sub_mult$BC_id, levels = paste0("BC_", 1:nrow(BC_list_num)))
cell_BC_sub_mult = cell_BC_sub_mult %>% group_by(cell_id) %>% arrange(BC_id, .by_group = T)
### 4%

cell_mult_BC = NULL
for(i in unique(cell_BC_sub_mult$cell_id)){
    x = as.vector(cell_BC_sub_mult$BC_id[cell_BC_sub_mult$cell_id == i])
    cell_mult_BC = rbind(cell_mult_BC, data.frame(cell_id = i,
                                                  group = paste(x, collapse = ",")))
}
cell_mult_BC %>% group_by(group) %>% tally() %>% arrange(desc(n))


#######################################
### Step-5: merging clones with MOI > 1

### BC_1 & BC_21; 
### BC_5 & BC_22;
### BC_3 & BC_24;
### BC_20 & BC_26

cell_mult_BC_sub = cell_mult_BC %>% filter(group %in% c("BC_1,BC_21", "BC_5,BC_22","BC_3,BC_24","BC_20,BC_26"))
cell_BC_sub_mult_sub = cell_BC_sub_mult %>% filter(cell_id %in% cell_mult_BC_sub$cell_id)

cell_BC_assigned = rbind(cell_BC_sub_uniq, cell_BC_sub_mult_sub)

BC = as.vector(cell_BC_assigned$BC)
BC_id = as.vector(cell_BC_assigned$BC_id)

x = paste0("BC_", c(1,21))
BC[cell_BC_assigned$BC_id %in% x] = paste(BC_list_num$BC[BC_list_num$BC_id %in% x], collapse = ",")
BC_id[cell_BC_assigned$BC_id %in% x] = paste(x, collapse = ",")

x = paste0("BC_", c(5,22))
BC[cell_BC_assigned$BC_id %in% x] = paste(BC_list_num$BC[BC_list_num$BC_id %in% x], collapse = ",")
BC_id[cell_BC_assigned$BC_id %in% x] = paste(x, collapse = ",")

x = paste0("BC_", c(3,24))
BC[cell_BC_assigned$BC_id %in% x] = paste(BC_list_num$BC[BC_list_num$BC_id %in% x], collapse = ",")
BC_id[cell_BC_assigned$BC_id %in% x] = paste(x, collapse = ",")

x = paste0("BC_", c(20,26))
BC[cell_BC_assigned$BC_id %in% x] = paste(BC_list_num$BC[BC_list_num$BC_id %in% x], collapse = ",")
BC_id[cell_BC_assigned$BC_id %in% x] = paste(x, collapse = ",")

cell_BC_assigned$BC = as.vector(BC)
cell_BC_assigned$BC_id = as.vector(BC_id)

cell_BC_assigned = cell_BC_assigned %>% select(BC, BC_id, cell_id) %>% unique() %>% as.data.frame()

clonotype = cell_BC_assigned %>% group_by(BC) %>% tally() %>% arrange(desc(n)) %>% filter(n >= 50)
clonotype$clonotype_id = paste0("clonotype_", 1:nrow(clonotype))

cell_BC_assigned = cell_BC_assigned %>%
    left_join(clonotype %>% select(BC, clonotype_id), by = "BC") %>%
    filter(!is.na(clonotype_id)) 


saveRDS(cell_BC_assigned, paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/cell_clonotype_rep1.rds"))
print(nrow(cell_BC_assigned)/sum(pd$experiment_id %in% paste0("monoclonal_EB_proof_E2E_rep", c("1_1","1_2"))))
### Finally, 83.3% of cells with transcriptome have been assigned to a distinct clonotype in rep1

print(length(unique(cell_BC_assigned$cell_id))/length(unique(cell_BC_sub$cell_id)))





#########################################################################
### Step-6: Plotting UMAP for cells before filtering out unassigned cells

dat_x = bc_matrix[,colnames(bc_matrix) %in% cell_BC_sub$cell_id]
### 26 x 16512 cells

pc_num = min(15, nrow(dat_x))

obj = CreateSeuratObject(dat_x)
obj = NormalizeData(obj, normalization.method = "RC", scale.factor = 10000)
obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = nrow(obj))
obj = ScaleData(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = pc_num, verbose = FALSE)
obj = RunUMAP(object = obj, reduction = "pca", dims = 1:pc_num, min.dist = 0.1, n.neighbors = 30, n.components = 2)

pd = data.frame(obj[[]])
pd$UMAP_1 = Embeddings(obj, reduction = "umap")[,1]
pd$UMAP_2 = Embeddings(obj, reduction = "umap")[,2]
pd$cell_id = rownames(pd)

cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/cell_clonotype_rep1.rds"))
pd = pd %>% left_join(cell_BC_assigned[,c("clonotype_id", "cell_id")]) ### 15870 are assigned

clonotype_color_map = c("#d374df",
                        "#64b948",
                        "#6f65d8",
                        "#b2b135",
                        "#a141a5",
                        "#5ebf82",
                        "#d94e94",
                        "#588131",
                        "#6a6fbb",
                        "#db923c",
                        "#5f9ed7",
                        "#cb542b",
                        "#3fc1bf",
                        "#d4445a",
                        "#39855f",
                        "#cd8bc5",
                        "#b7aa65",
                        "#9d496c",
                        "#876a2b",
                        "#cc7667")
names(clonotype_color_map) = paste0("clonotype_", 1:20)

pd_x = pd %>% filter(!is.na(clonotype_id)) %>%
    group_by(clonotype_id) %>% summarize(UMAP_1_mean = mean(UMAP_1), UMAP_2_mean = mean(UMAP_2))
pd_x$clonotype_id = gsub("clonotype_", "", pd_x$clonotype_id)

p = ggplot() +
        geom_point(data = pd %>% filter(is.na(clonotype_id)), aes(x = UMAP_1, y = UMAP_2), color = "grey70", size=0.2) +
        geom_point(data = pd %>% filter(!is.na(clonotype_id)), aes(x = UMAP_1, y = UMAP_2, color = clonotype_id), size=0.2) +
        ggrepel::geom_text_repel(data = pd_x, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = clonotype_id), color = "black", size = 3, family = "Arial") +
        theme_void() +
        scale_color_manual(values=clonotype_color_map) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/UMAP_clonotype_before_filtering_rep1.png"), p,
               dpi = 300,
               height  = 5, 
               width = 5)




################################################################
### Step-7: Identifying enriched gRNAs for individual clonotypes

gRNA_matrix = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/gRNA_matrix_rep1.rds"))
cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/cell_clonotype_rep1.rds"))
gRNA_matrix = gRNA_matrix[,colnames(gRNA_matrix) %in% as.vector(cell_BC_assigned$cell_id)]
gRNA_matrix = gRNA_matrix[,as.vector(cell_BC_assigned$cell_id)]

clonotype_list = unique(cell_BC_assigned$clonotype_id)
cell_gRNA = NULL
for(i in clonotype_list){
    gRNA_matrix_i = gRNA_matrix[,cell_BC_assigned$clonotype_id == i]
    gRNA_matrix_i = gRNA_matrix_i[,colSums(gRNA_matrix_i) >= 10]
    gRNA_matrix_i_norm = t(t(gRNA_matrix_i)/colSums(gRNA_matrix_i))
    gRNA_matrix_i_norm_rowMean = sort(rowMeans(gRNA_matrix_i_norm), decreasing=T)
    
    cell_gRNA = rbind(cell_gRNA, data.frame(clonotype_id = i,
                                            gRNA = names(gRNA_matrix_i_norm_rowMean),
                                            pct = 100*as.vector(gRNA_matrix_i_norm_rowMean),
                                            order = 1:length(gRNA_matrix_i_norm_rowMean)))
}

saveRDS(cell_gRNA, paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/cell_clonotype_gRNA_rep1.rds"))




####################################################################################
### Step-8: Making a plot to present the distribution of UMIs of gRNAs per clonotype

cell_gRNA = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/cell_clonotype_gRNA_rep1.rds"))

df = NULL
clonotype_list = unique(cell_gRNA$clonotype_id)
for(i in clonotype_list){
    df_i = subset(cell_gRNA, clonotype_id == i)
    df_i$pct_cum = cumsum(df_i$pct)
    df = rbind(df, df_i)
}

ggplot(df, aes(x = order, y = pct_cum, color = clonotype_id)) + 
    geom_point() + 
    geom_line(linewidth = 1) +
    labs(title = "", x = "gRNAs", y = "cumulative % of UMIs") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    scale_color_manual(values=clonotype_color_map) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/cell_clonotype_gRNA_rep1_line_plot.pdf"),
           dpi = 300,
           height  = 5, 
           width = 5)


cell_gRNA_sub = cell_gRNA %>% filter(pct >= 10) %>% arrange(clonotype_id)


################################################################################
#### Step-9: validating using the "look-up" table; of note, this "look-up" table is based on sequencing result of the library, rather than experimental design

library(stringr)
white_list = read.table("NTC_scaled_sgRNA_BC_final_list_20210510.txt", header=T, as.is=T)
white_list$old_sgRNA = paste0("G", white_list$sgRNA)
white_list$sgRNA = paste0("G", substr(white_list$sgRNA, 1, nchar(white_list$sgRNA)-1))

cell_gRNA = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/cell_clonotype_gRNA_rep1.rds"))
cell_gRNA = cell_gRNA[cell_gRNA$pct >= 10,]
cell_gRNA$gRNA = unlist(lapply(as.vector(cell_gRNA$gRNA), function(x) strsplit(x,"[_]")[[1]][2])) 

cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/cell_clonotype_rep1.rds"))
cell_BC_assigned = unique(cell_BC_assigned[,c("BC", "clonotype_id")])
df = NULL
for(i in 1:nrow(cell_BC_assigned)){
    if(grepl(',', cell_BC_assigned$BC[i])){
        x1 = unlist(lapply(as.vector(cell_BC_assigned$BC[i]), function(x) strsplit(x,"[,]")[[1]][1])) 
        df = rbind(df, data.frame(BC = x1, clonotype_id = cell_BC_assigned$clonotype_id[i]))
        x2 = unlist(lapply(as.vector(cell_BC_assigned$BC[i]), function(x) strsplit(x,"[,]")[[1]][2])) 
        df = rbind(df, data.frame(BC = x2, clonotype_id = cell_BC_assigned$clonotype_id[i]))
    } else {
        df = rbind(df, data.frame(BC = cell_BC_assigned$BC[i], clonotype_id = cell_BC_assigned$clonotype_id[i]))
    }
}

reverse_complement <- function(dna_sequences) {
    complement <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C", "N" = "N")
    reverse_comp_single <- function(sequence) {
        paste(rev(complement[strsplit(sequence, split = "")[[1]]]), collapse = "")
    }
    sapply(dna_sequences, reverse_comp_single)
}
df$BC_rev = as.vector(reverse_complement(df$BC))
df = df %>% left_join(cell_gRNA[,c(1:2)], by = "clonotype_id", relationship = "many-to-many")

res = df %>% select(sgRNA = gRNA, BC = BC_rev) %>%
    left_join(white_list, by = c("sgRNA", "BC"))
sum(!is.na(res$sgRNA_BC_read_count))
### n = 13 (out of 31)

res_2 = df %>% select(old_sgRNA = gRNA, BC = BC_rev) %>%
    left_join(white_list, by = c("old_sgRNA", "BC"))
sum(!is.na(res_2$sgRNA_BC_read_count))
### n = 4


###########################################################
### Step-10: Making separated UMAP on the top 10 clonotypes


cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/cell_clonotype_rep1.rds"))
pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_proof/obj_processed_pd.rds"))

for(i in paste0("clonotype_", 1:10)){
    try(ggplot() +
            geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "grey80", size=0.1) +
            geom_point(data = pd %>% filter(cell_id %in% as.vector(cell_BC_assigned$cell_id[cell_BC_assigned$clonotype_id == i])), aes(x = UMAP_1, y = UMAP_2), size=0.5, color = "red") +
            theme_void() +
            theme(legend.position="none") +
            theme(plot.title = element_text(hjust = 0.5)) +
            ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/UMAP_rep1_",  i, ".png"),
                   dpi = 300,
                   height  = 5, 
                   width = 5), silent = TRUE)
}

df = cell_BC_assigned %>% 
    left_join(pd %>% select(cell_id, celltype), by = "cell_id") %>%
    group_by(clonotype_id, celltype) %>%
    tally()
df$clonotype_id = factor(df$clonotype_id, levels = paste0("clonotype_", 1:20))
df$celltype = factor(df$celltype, levels = c("Primordial germ cells",
                                             "Epiblast",
                                             "Primitive streak",
                                             "Neuroectoderm",
                                             "Surface ectoderm",
                                             "Notochord",
                                             "Definitive endoderm",
                                             "ExE visceral endoderm",
                                             "Paraxial mesoderm",
                                             "Nascent mesoderm",
                                             "Lateral plate mesoderm",
                                             "Cardiomyocytes",
                                             "Endothelial cells",
                                             "Blood progenitors",
                                             "Erythroid cells",
                                             "Early neurons"))

# Stacked + percent across clonotypes
p = df %>%
    ggplot(aes(fill=celltype, y=n, x=clonotype_id)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=EB_celltype_color_code) +
    labs(x="", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/celltype_compositions_rep1.pdf"),
           dpi = 300,
           height  = 5, 
           width = 5)

clonotype_cell_num = cell_BC_assigned %>% group_by(clonotype_id) %>% tally() %>% arrange(desc(n))
clonotype_cell_num$log2_cell_num = log2(clonotype_cell_num$n)
clonotype_cell_num$clonotype_id = factor(clonotype_cell_num$clonotype_id, levels = as.vector(clonotype_cell_num$clonotype_id))

p = ggplot(data=clonotype_cell_num, aes(x=clonotype_id, y=n)) +
    geom_bar(stat="identity") + theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_fill_viridis() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    #scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = "cell_count")) +
    #guides(y = "none") +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_proof/cell_to_mEB/clonotype_cell_num_rep1.pdf"),
           height  = 2, 
           width = 5)










