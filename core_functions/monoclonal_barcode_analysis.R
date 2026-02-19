
#############################################################################
### Assiging cells to each mEBs based on the barcode, on "arrayed" experiment
### contact: cxqiu@uw.edu

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))

run_list = paste0("clonal_EBTF_", c(1:12))

#############################
### Step-1: reading BC matrix 

bc_matrix = NULL

for(run_id in run_list){
    print(run_id)
    if(run_id %in% paste0("clonal_EBTF_", c(1,2,3,4))){
        bc_matrix_table = read.table(paste0(work_path, "/data/monoclonal_EB_TF_screen/clonal_EBTF_1_thru_4_get_bc.txt"), header=T, as.is=T)
    } else if (run_id %in% paste0("clonal_EBTF_", c(5,6,7,8))){
        bc_matrix_table = read.table(paste0(work_path, "/data/monoclonal_EB_TF_screen/clonal_EBTF_5_thru_8_get_bc.txt"), header=T, as.is=T)
    } else {
        bc_matrix_table = read.table(paste0(work_path, "/data/monoclonal_EB_TF_screen/clonal_EBTF_9_thru_12_get_bc.txt"), header=T, as.is=T)
    }
    bc_matrix_table$cell_id = paste0("monoclonal_EB_TF_screen_", run_id, "_",  unlist(lapply(as.vector(bc_matrix_table$cBC), function(x) strsplit(x,"[-]")[[1]][1])))
    bc_matrix_table = subset(bc_matrix_table, cell_id %in% as.vector(pd$cell_id))
    
    bc_matrix = rbind(bc_matrix, bc_matrix_table)
}

########################################################################
### Step-2: performing barcode correction based on hamming distance == 1

tmp_rep = bc_matrix %>% group_by(mBC) %>%
    summarize(read_count_sum = sum(n_reads_filtered), umi_count_sum = sum(filtered_corrected_UMIs)) %>%
    filter(umi_count_sum >= 10, read_count_sum >= 10) %>% arrange(desc(umi_count_sum))
tmp_rep$mBC_id = paste0("BC_", 1:nrow(tmp_rep))

hamming_distance <- function(string1, string2) {
    # Split the strings into individual characters
    chars1 <- strsplit(string1, NULL)[[1]]
    chars2 <- strsplit(string2, NULL)[[1]]
    
    # Compare the characters and sum the number of differences
    distance <- sum(chars1 != chars2)
    
    return(distance)
}

tmp_rep_sequence = list()
for(i in 1:nrow(tmp_rep)){
    tmp_rep_sequence[[i]] = strsplit(tmp_rep$mBC[i], NULL)[[1]]
}

start_time = Sys.time()
ham_dis_bc = NULL
for(i in 1:(nrow(tmp_rep)-1)){
    for(j in (i+1):nrow(tmp_rep)){
        if(sum(tmp_rep_sequence[[i]] != tmp_rep_sequence[[j]]) == 1){
            ham_dis_bc = c(ham_dis_bc, c(i, j))
        }
    }
}
print(Sys.time() - start_time)

edges = paste0("BC_", ham_dis_bc)
g = graph(edges, directed = FALSE)
clusters_info = clusters(g)
components = split(names(clusters_info$membership), clusters_info$membership)

bc_correction_list = NULL
for(i in 1:length(components)){
    print(i)
    y = paste0("BC_", 1:nrow(tmp_rep))[paste0("BC_", 1:nrow(tmp_rep)) %in% components[[i]]]
    for(j in y){
        bc_correction_list = rbind(bc_correction_list, data.frame(mBC_id = j, mBC_corrected = y[1]))
    }
}

tmp_rep_1 = tmp_rep %>% filter(mBC_id %in% bc_correction_list$mBC_id) %>%
    left_join(bc_correction_list, by = "mBC_id")
tmp_rep_2 = tmp_rep %>% filter(!mBC_id %in% bc_correction_list$mBC_id) %>%
    mutate(mBC_corrected = mBC_id)
tmp_rep_x = rbind(tmp_rep_1, tmp_rep_2)
tmp_rep_x = tmp_rep_x %>% left_join(tmp_rep %>% select(mBC_corrected = mBC_id, mBC_correct_seq = mBC), by = "mBC_corrected") %>%
    as.data.frame()

tmp_rep_summary = tmp_rep_x %>% group_by(mBC_correct_seq) %>% 
    summarize(read_count_sum_corrected = sum(read_count_sum), umi_count_sum_corrected = sum(umi_count_sum)) %>%
    mutate(log2_read_count = log2(read_count_sum_corrected), log2_umi_count = log2(umi_count_sum_corrected))

try(ggplot() +
        geom_point(data = subset(tmp_rep_summary, log2_read_count < 15), aes(x = log2_read_count, y = log2_umi_count), size = 0.5) + 
        geom_point(data = subset(tmp_rep_summary, log2_read_count >= 15), aes(x = log2_read_count, y = log2_umi_count), color = "red", size = 1) +
        geom_vline(xintercept = 15) +
        labs(x="Log2 read count per BC", y="Log2 UMI count per BC", title="") +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/BC_2D_hist_rep1.pdf"),
               height  = 3, 
               width = 5), silent = TRUE)


### After checking the 2D-hist, I decided to use log2_read_count >= 15 as cutoff to filter BCs
### n = 331 BCs are retained
tmp_include = subset(tmp_rep_summary, log2_read_count >= 15)

### filtering BCs and gRNAs
bc_matrix_list = list()
gRNA_matrix_list = list()

for(run_id in run_list){    
    print(run_id)
    pd_filter = subset(pd, experiment_id == paste0("monoclonal_EB_TF_screen_", run_id))
    
    count_matrix_all = Read10X(paste0(work_path, "/data/monoclonal_EB_TF_screen/", run_id, "/outs/raw_feature_bc_matrix"), 
                               gene.column = 1,
                               strip.suffix = T)
    gRNA_matrix = count_matrix_all[["CRISPR Guide Capture"]]
    colnames(gRNA_matrix) = paste0("monoclonal_EB_TF_screen_", run_id, "_", colnames(gRNA_matrix))
    gRNA_matrix = gRNA_matrix[,colnames(gRNA_matrix) %in% as.vector(pd_filter$cell_id)]
    gRNA_matrix = gRNA_matrix[rowMeans(gRNA_matrix > 0) >= 0.01,]
    
    gRNA_matrix_list[[run_id]] = gRNA_matrix
    
    if(run_id %in% paste0("clonal_EBTF_", c(1,2,3,4))){
        bc_matrix_table = read.table(paste0(work_path, "/data/monoclonal_EB_TF_screen/clonal_EBTF_1_thru_4_get_bc.txt"), header=T, as.is=T)
    } else if (run_id %in% paste0("clonal_EBTF_", c(5,6,7,8))){
        bc_matrix_table = read.table(paste0(work_path, "/data/monoclonal_EB_TF_screen/clonal_EBTF_5_thru_8_get_bc.txt"), header=T, as.is=T)
    } else {
        bc_matrix_table = read.table(paste0(work_path, "/data/monoclonal_EB_TF_screen/clonal_EBTF_9_thru_12_get_bc.txt"), header=T, as.is=T)
    }
    bc_matrix_table$cell_id = paste0("monoclonal_EB_TF_screen_", run_id, "_",  unlist(lapply(as.vector(bc_matrix_table$cBC), function(x) strsplit(x,"[-]")[[1]][1])))
    bc_matrix_table = subset(bc_matrix_table, cell_id %in% as.vector(pd_filter$cell_id) & mBC %in% as.vector(tmp_rep_x$mBC))
    
    bc_matrix_table_x = bc_matrix_table %>% left_join(tmp_rep_x[,c("mBC","mBC_correct_seq")], by = "mBC") %>% as.data.frame()
    bc_matrix_table_x = bc_matrix_table_x %>% group_by(mBC_correct_seq, cell_id) %>% summarize(filtered_corrected_UMIs_x = sum(filtered_corrected_UMIs))
    
    bc_matrix = bc_matrix_table_x[,c("mBC_correct_seq", "cell_id", "filtered_corrected_UMIs_x")]
    bc_matrix = dcast(bc_matrix, mBC_correct_seq ~ cell_id, fill = 0)
    rownames(bc_matrix) = bc_matrix[,1]
    bc_matrix = bc_matrix[,-1]
    bc_matrix = as(bc_matrix, "sparseMatrix")
    bc_matrix = bc_matrix[rownames(bc_matrix) %in% as.vector(tmp_include$mBC_correct_seq),]
    
    bc_matrix_list[[run_id]] = bc_matrix
}

### double-checking if feature names are matched between replicates
### n = 331 BCs, n = 30 gRNAs
print(sum(rownames(bc_matrix_list[[1]]) %in% rownames(bc_matrix_list[[2]])))
print(sum(rownames(gRNA_matrix_list[[1]]) %in% rownames(gRNA_matrix_list[[2]])))

bc_matrix = NULL
gRNA_matrix = NULL

bc_matrix_row = unique(tmp_include$mBC_correct_seq)
gRNA_matrix_row = rownames(gRNA_matrix_list[[1]])

for(run_id in run_list){
    print(run_id)
    tmp = bc_matrix_list[[run_id]]
    if(nrow(tmp) < length(bc_matrix_row)){
        feature_not_include = bc_matrix_row[!bc_matrix_row %in% rownames(tmp)]
        tmp_add = matrix(0, length(feature_not_include), ncol(tmp))
        rownames(tmp_add) = feature_not_include
        colnames(tmp_add) = colnames(tmp)
        tmp = rbind(tmp, tmp_add)
    }
    
    if(is.null(bc_matrix)){
        bc_matrix = tmp
    } else {
        bc_matrix = cbind(bc_matrix, tmp[rownames(bc_matrix),])
    }
    
    tmp = gRNA_matrix_list[[run_id]]
    if(nrow(tmp) < length(gRNA_matrix_row)){
        feature_not_include = gRNA_matrix_row[!gRNA_matrix_row %in% rownames(tmp)]
        tmp_add = matrix(0, length(feature_not_include), ncol(tmp))
        rownames(tmp_add) = feature_not_include
        colnames(tmp_add) = colnames(tmp)
        tmp = rbind(tmp, tmp_add)
    }
    
    if(is.null(gRNA_matrix)){
        gRNA_matrix = tmp
    } else {
        gRNA_matrix = cbind(gRNA_matrix, tmp[rownames(gRNA_matrix),])
    }
}

saveRDS(bc_matrix, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/BC_matrix_rep1.rds"))
saveRDS(gRNA_matrix, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/gRNA_matrix_rep1.rds"))



#########################################################################
### Step-3: Setting cutoffs to assign cells to clonotype based on the BCs

cell_BC = 100*t(t(bc_matrix)/colSums(bc_matrix))
cell_BC = melt(as.matrix(cell_BC))
names(cell_BC) = c("BC", "cell_id", "pct")
cell_BC_x = melt(as.matrix(bc_matrix))
cell_BC$UMI_count = as.vector(cell_BC_x$value)
cell_BC = cell_BC %>% group_by(cell_id) %>% arrange(desc(pct), .by_group = TRUE)
cell_BC$pct_order = rep(c(1:nrow(bc_matrix)), times = ncol(bc_matrix))

ggplot(data = subset(cell_BC, UMI_count > 0), aes(log2(UMI_count))) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = log2(10)) + theme_classic(base_size = 10) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/hist_umi_per_BC_per_cell_rep1.pdf"),
           height  = 4, 
           width = 4)

ggplot(data = subset(cell_BC, UMI_count > 0), aes(pct)) + geom_histogram(bins = 30) + 
    geom_vline(xintercept = 30) + theme_classic(base_size = 10) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/hist_pct_umi_per_BC_per_cell_rep1.pdf"),
           height  = 4, 
           width = 4)


### pct >= 30 and UMI_count >= 10
cell_BC_sub = subset(cell_BC, pct >= 30 & UMI_count >= 10)
print(length(unique(cell_BC_sub$cell_id))/length(unique(cell_BC$cell_id)))
print(sum(table(cell_BC_sub$cell_id) == 1)/length(unique(cell_BC$cell_id)))
### 85.6% cells are retained, and 73.9% cells have unique BC




##############################################################
### Step-4: calculating Peasrson correlation on the UMI matrix

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
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/BC_cell_num.pdf"),
           height  = 5, 
           width = 5)


### calculating Peasrson correlation on the UMI matrix
bc_matrix_sub = bc_matrix[,colnames(bc_matrix) %in% cell_BC_sub$cell_id]
bc_matrix_sub_norm = t(bc_matrix_sub)/colSums(bc_matrix_sub)
bc_matrix_sub_scale = scale(bc_matrix_sub_norm)
bc_matrix_sub_scale = bc_matrix_sub_scale[,BC_list]

dat = cor(bc_matrix_sub_scale)
dat = data.frame(A = rep(1:nrow(dat), times = nrow(dat)),
                 B = rep(1:nrow(dat), each = nrow(dat)),
                 corr = c(dat))
dat = subset(dat, A >= B)

p = ggplot() + 
    geom_point(data = subset(dat, A == B), aes(A, B), shape=22, color="grey80", fill="white") + 
    geom_point(data = subset(dat, A != B), aes(A, B, fill = corr), shape=22, color="grey80") + 
    labs(x="", y="", title="") +
    theme_classic(base_size = 10) +
    #theme(legend.position="none") +
    scale_fill_viridis() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/heatmap_UMI_correlation_BC.pdf"),
           height  = 15, 
           width = 18)

### identifying clones with MOI > 1
dat_x = dat %>% filter(A != B, corr > 0.1) %>% arrange(desc(corr)) 


###########################################
### Step-5: identifying distinct clonotypes

cell_BC_sub_uniq_x = cell_BC_sub %>% group_by(cell_id) %>% tally() %>% filter(n == 1)

cell_BC_sub_uniq = cell_BC_sub %>% filter(cell_id %in% as.vector(cell_BC_sub_uniq_x$cell_id)) %>%
    left_join(BC_list_num[,c("BC","BC_id")], by = "BC")

cell_BC_sub_mult = cell_BC_sub %>% filter(!cell_id %in% as.vector(cell_BC_sub_uniq_x$cell_id)) %>%
    left_join(BC_list_num[,c("BC","BC_id")], by = "BC")
cell_BC_sub_mult$BC_id = factor(cell_BC_sub_mult$BC_id, levels = paste0("BC_", 1:nrow(BC_list_num)))
cell_BC_sub_mult = cell_BC_sub_mult %>% group_by(cell_id) %>% arrange(BC_id, .by_group = T)

cell_mult_BC = NULL
for(i in unique(cell_BC_sub_mult$cell_id)){
    x = as.vector(cell_BC_sub_mult$BC_id[cell_BC_sub_mult$cell_id == i])
    cell_mult_BC = rbind(cell_mult_BC, data.frame(cell_id = i,
                                                  group = paste(x, collapse = ",")))
}
cell_mult_BC %>% group_by(group) %>% tally() %>% arrange(desc(n))


#######################################
### Step-6: merging clones with MOI > 1

library(igraph)
edges = paste0("BC_", c(t(as.matrix(dat_x[,c(1,2)]))))
g = graph(edges, directed = FALSE)
clusters_info = clusters(g)
components = split(names(clusters_info$membership), clusters_info$membership)

# Function to generate subsets of size greater than 1
generate_subsets <- function(x) {
    # Get all subsets of size greater than 1
    subsets <- lapply(2:length(x), function(i) combn(x, i, simplify = FALSE))
    # Flatten the list of lists
    subsets <- unlist(subsets, recursive = FALSE)
    return(subsets)
}

x_list = NULL
for(i in 1:length(components)){
    y = generate_subsets(BC_list_num$BC_id[BC_list_num$BC_id %in% components[[i]]])
    for(j in 1:length(y)){
        x_list = rbind(x_list, 
                       data.frame(A = paste(y[[j]], collapse = ","),
                                  B = paste(BC_list_num$BC_id[BC_list_num$BC_id %in% components[[i]]], collapse = ",")))
    }
}

cell_mult_BC_sub = cell_mult_BC %>% filter(group %in% as.vector(x_list$A))
cell_BC_sub_mult_sub = cell_BC_sub_mult %>% filter(cell_id %in% cell_mult_BC_sub$cell_id)

cell_BC_assigned = rbind(cell_BC_sub_uniq, cell_BC_sub_mult_sub)

BC = as.vector(cell_BC_assigned$BC)
BC_id = as.vector(cell_BC_assigned$BC_id)

for(i in 1:length(components)){
    y = BC_list_num$BC_id[BC_list_num$BC_id %in% components[[i]]]
    BC[cell_BC_assigned$BC_id %in% y] = paste(BC_list_num$BC[BC_list_num$BC_id %in% components[[i]]], collapse = ",")
    BC_id[cell_BC_assigned$BC_id %in% y] = paste(y, collapse = ",")
}

cell_BC_assigned$BC = as.vector(BC)
cell_BC_assigned$BC_id = as.vector(BC_id)

cell_BC_assigned = cell_BC_assigned %>% select(BC, BC_id, cell_id) %>% unique() %>% as.data.frame()

clonotype = cell_BC_assigned %>% group_by(BC) %>% tally() %>% arrange(desc(n)) %>% filter(n >= 50)
clonotype$clonotype_id = paste0("clonotype_", 1:nrow(clonotype))

cell_BC_assigned = cell_BC_assigned %>%
    left_join(clonotype %>% select(BC, clonotype_id), by = "BC") %>%
    filter(!is.na(clonotype_id)) 



saveRDS(cell_BC_assigned, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep1.rds"))
print(nrow(cell_BC_assigned)/sum(pd$experiment_id %in% paste0("monoclonal_EB_TF_screen_clonal_EBTF_", c(1:12))))
### Finally, 79.1% of cells with transcriptome have been assigned to a distinct clonotype in rep1


clonotype$clonotype_id = factor(clonotype$clonotype_id, levels = paste0("clonotype_", 1:nrow(clonotype)))
p = ggplot(data=clonotype, aes(x=clonotype_id, y=n)) +
    geom_bar(stat="identity") + theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_fill_viridis() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    #scale_y_continuous(trans='log2') +
    #scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = "cell_count")) +
    #guides(y = "none") +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/clonotype_cell_num.pdf"),
           height  = 2.5, 
           width = 5)



#########################################################################
### Step-7: Plotting UMAP for cells before filtering out unassigned cells

dat_x = bc_matrix[,colnames(bc_matrix) %in% cell_BC_sub$cell_id]
### 331 x 58,232

pc_num = min(30, nrow(dat_x))

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

cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep1.rds"))
pd = pd %>% left_join(cell_BC_assigned[,c("clonotype_id", "cell_id")]) ### 55063 are assigned

try(ggplot() +
        geom_point(data = pd %>% filter(is.na(clonotype_id)), aes(x = UMAP_1, y = UMAP_2), color = "grey70", size=0.2) +
        geom_point(data = pd %>% filter(!is.na(clonotype_id)), aes(x = UMAP_1, y = UMAP_2, color = clonotype_id), size=0.2) +
        # ggrepel::geom_text_repel(data = cell_BC_assigned %>% group_by(clonotype_id) %>% sample_n(1), aes(x = UMAP_1, y = UMAP_2, label = clonotype), color = "black", size = 3, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/UMAP_clonotype_before_filtering_rep1_nolabel.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

################################################################
### Step-8: Identifying enriched gRNAs for individual clonotypes

gRNA_matrix = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/gRNA_matrix_rep1.rds"))
cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep1.rds"))
gRNA_matrix = gRNA_matrix[,colnames(gRNA_matrix) %in% as.vector(cell_BC_assigned$cell_id)]
gRNA_matrix = gRNA_matrix[,as.vector(cell_BC_assigned$cell_id)]

clonotype_list = unique(cell_BC_assigned$clonotype_id)
cell_gRNA = NULL
for(i in clonotype_list){
    gRNA_matrix_i = gRNA_matrix[,cell_BC_assigned$clonotype_id == i]
    if(sum(colSums(gRNA_matrix_i) >= 10) >= 10){
        gRNA_matrix_i = gRNA_matrix_i[,colSums(gRNA_matrix_i) >= 10]
        gRNA_matrix_i_norm = t(t(gRNA_matrix_i)/colSums(gRNA_matrix_i))
        gRNA_matrix_i_norm_rowMean = sort(rowMeans(gRNA_matrix_i_norm), decreasing=T)
        
        cell_gRNA = rbind(cell_gRNA, data.frame(clonotype_id = i,
                                                gRNA = names(gRNA_matrix_i_norm_rowMean),
                                                pct = 100*as.vector(gRNA_matrix_i_norm_rowMean),
                                                order = 1:length(gRNA_matrix_i_norm_rowMean)))
    }
}

saveRDS(cell_gRNA, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep1.rds"))


####################################################################################
### Step-9: Making a plot to present the distribution of UMIs of gRNAs per clonotype

cell_gRNA = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep1.rds"))

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
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep1_line_plot.pdf"),
           dpi = 300,
           height  = 5, 
           width = 5)


cell_gRNA_sub = cell_gRNA %>% filter(pct >= 50) %>% arrange(clonotype_id)
cell_gRNA_sub$target = unlist(lapply(as.vector(cell_gRNA_sub$gRNA), function(x) strsplit(x,"[_]")[[1]][1])) 



############################################################################################################
### Step-10: Checking the pair of oBC and gRNA identified from the data is consistent with the look up table

### gRNA_BC look up table
gRNA_BC = read.table("/net/shendure/vol8/projects/cxqiu/nobackup/work/sam_tf/processing/monoclonal_EB_TF_screen_sci/gRNA_BC_look_up_table.txt", as.is=T)
colnames(gRNA_BC) = c("gRNA", "BC")

### gRNA assignments for individual clonotype
cell_gRNA = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep1.rds"))

### oBC assignments for individual clonotype
bc_matrix = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/BC_matrix_rep1.rds"))
cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep1.rds"))
bc_matrix = bc_matrix[,colnames(bc_matrix) %in% as.vector(cell_BC_assigned$cell_id)]
bc_matrix = bc_matrix[,as.vector(cell_BC_assigned$cell_id)]

clonotype_list = unique(cell_BC_assigned$clonotype_id)
cell_bc = NULL
for(i in clonotype_list){
    bc_matrix_i = bc_matrix[,cell_BC_assigned$clonotype_id == i]
    if(sum(colSums(bc_matrix_i) >= 10) >= 10){
        bc_matrix_i = bc_matrix_i[,colSums(bc_matrix_i) >= 10]
        bc_matrix_i_norm = t(t(bc_matrix_i)/colSums(bc_matrix_i))
        bc_matrix_i_norm_rowMean = sort(rowMeans(bc_matrix_i_norm), decreasing=T)
        
        cell_bc = rbind(cell_bc, data.frame(clonotype_id = i,
                                            BC = names(bc_matrix_i_norm_rowMean),
                                            pct = 100*as.vector(bc_matrix_i_norm_rowMean),
                                            order = 1:length(bc_matrix_i_norm_rowMean)))
    }
}

### For individual clonotype, plot the top abundant oBC vs. its matched gRNA from the look up table
dat = cell_bc %>% filter(order == 1) %>% select(clonotype_id, org_BC = BC, pct) %>%
    mutate(BC = stringr::str_sub(org_BC, 2, 9)) %>%
    left_join(gRNA_BC, by = "BC") 
dat$gRNA = gsub('[.]','_',dat$gRNA)
dat = dat %>% left_join(cell_gRNA %>% select(clonotype_id, gRNA, gRNA_pct = pct), by = c("gRNA","clonotype_id"))

dat = dat[!is.na(dat$gRNA_pct),]
sum(dat$pct >= 85 & dat$gRNA_pct >= 85)/nrow(dat) ### 80.7%

### separately color uniq assigned clone or double assigned
### out of 202 clonotypes, 187 are assigned with uique gRNAs, 15 are assigned with double gRNAs
cell_gRNA_uniq = cell_gRNA %>% filter(pct >= 70) %>% arrange(clonotype_id)
dat$uniq_mult = if_else(dat$clonotype_id %in% as.vector(cell_gRNA_uniq$clonotype_id), "unique", "double")

p = dat %>% 
    ggplot(aes(pct, gRNA_pct, color = uniq_mult)) + geom_point(size=0.5) + 
    labs(x="% of UMIs from the top oBC", y="% of UMIs from the paired gRNA", title="") +
    xlim(0, 100) + ylim(0, 100) +
    geom_vline(xintercept = 85, linetype = "dotted") +
    geom_hline(yintercept = 85, linetype = "dotted") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_color_manual(values=c("unique" = "#FF0000", "double" = "#0000FF")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
ggsave("~/share/UMIs_of_oBC_gRNA_pairs.pdf", p, width=3, height=3)




