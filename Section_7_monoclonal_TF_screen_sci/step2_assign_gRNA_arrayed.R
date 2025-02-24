
##########################################################
### Assigning sgRNA to monoclonal EBs (arrayed experiment)
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
library(igraph)

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen_sci/obj_processed_pd.rds"))
rownames(pd) = NULL
pd$cell_id = paste0(pd$PCR_p5_sequence, "_", pd$RT_barcode_well, "_", pd$LIG_barcode_well)

### BC UMI data
dat_BC_1 = read.table(paste0(data_path, "/ebBC.mEB.ebTF.SR.filtered.handpool.1to10-L001_sci_mBC_UMI_corrected_counts_20240726.txt.gz"), header=T, as.is=T)
dat_BC_1$RT_group = "GEx_mEB_ebTF_SR_filtered_handpool_1to10"
dat_BC_2 = read.table(paste0(data_path, "/ebBC.mEB.ebTF.SR.filtered.pool11-L001_sci_mBC_UMI_corrected_counts_20240726.txt.gz"), header=T, as.is=T)
dat_BC_2$RT_group = "GEx_mEB_ebTF_SR_filtered_pool11"
dat_BC_3 = read.table(paste0(data_path, "/ebBC.mEB.ebTF.SR.gravity.handpool.1to10-L001_sci_mBC_UMI_corrected_counts_20240726.txt.gz"), header=T, as.is=T)
dat_BC_3$RT_group = "GEx_mEB_ebTF_SR_gravity_handpool_1to10"
dat_BC_4 = read.table(paste0(data_path, "/ebBC.mEB.ebTF.SR.gravity.pool11-L001_sci_mBC_UMI_corrected_counts_20240726.txt.gz"), header=T, as.is=T)
dat_BC_4$RT_group = "GEx_mEB_ebTF_SR_gravity_pool11"

dat_BC = rbind(dat_BC_1, dat_BC_3)
tmp_1 = unlist(lapply(as.vector(dat_BC$cBC), function(x) strsplit(x,"[_]")[[1]][1])) 
tmp_2 = unlist(lapply(as.vector(dat_BC$cBC), function(x) strsplit(x,"[_]")[[1]][3])) 
tmp_3 = unlist(lapply(as.vector(dat_BC$cBC), function(x) strsplit(x,"[_]")[[1]][4])) 
dat_BC$cell_id = paste0(tmp_1, "_", tmp_2, "_", tmp_3)

### merging cells (JB designed replicates on P7 index, so merging them)
dat_BC = dat_BC %>% rename(mBC = BC_oi) %>% group_by(cell_id, mBC) %>% 
    summarize(umi_count = sum(filtered_corrected_UMIs),
              read_count = sum(n_reads_filtered))

### filtering cells with transcriptome
dat_BC = dat_BC %>% left_join(pd %>% select(cell_id, sample), by = "cell_id") %>%
    filter(!is.na(sample)) %>%
    as.data.frame()

### barcode correction
tmp_rep = dat_BC %>% group_by(mBC) %>%
    summarize(read_count_sum = sum(read_count), umi_count_sum = sum(umi_count)) %>%
    filter(umi_count_sum >= 10, read_count_sum >= 10) %>% arrange(desc(umi_count_sum))
tmp_rep$mBC_id = paste0("BC_", 1:nrow(tmp_rep))

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
        geom_vline(xintercept = 10) +
        labs(x="Log2 read count per BC", y="Log2 UMI count per BC", title="") +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen_sci/cell_to_mEB/BC_2D_hist_rep1.pdf"),
               height  = 3, 
               width = 5), silent = TRUE)

### After checking the 2D-hist, I decided to use log2_read_count >= 10 as cutoff to filter BCs
### n = 350 BCs are retained
tmp_include = subset(tmp_rep_summary, log2_read_count >= 10)

dat_BC_sub = dat_BC %>% filter(mBC %in% tmp_rep_x$mBC) %>%
    left_join(tmp_rep_x[,c("mBC", "mBC_correct_seq")], by = "mBC") %>%
    select(sample, mBC_correct_seq, umi_count) %>%
    rename(cell_id = sample) %>%
    group_by(cell_id, mBC_correct_seq) %>%
    summarize(umi_count_sum = sum(umi_count))

bc_matrix = dat_BC_sub[,c("mBC_correct_seq", "cell_id", "umi_count_sum")]
bc_matrix = dcast(bc_matrix, mBC_correct_seq ~ cell_id, fill = 0)
rownames(bc_matrix) = bc_matrix[,1]
bc_matrix = bc_matrix[,-1]
bc_matrix = as(bc_matrix, "sparseMatrix")
bc_matrix = bc_matrix[rownames(bc_matrix) %in% as.vector(tmp_include$mBC_correct_seq),]
bc_matrix = bc_matrix[,colSums(bc_matrix) >= 10]

saveRDS(bc_matrix, paste0(work_path, "/analysis/monoclonal_EB_TF_screen_sci/cell_to_mEB/bc_matrix_rep1.rds"))



### Setting cutoffs to assign cells to clonotype based on the BCs

cell_BC = 100*t(t(bc_matrix)/colSums(bc_matrix))
cell_BC = melt(as.matrix(cell_BC))
names(cell_BC) = c("BC", "cell_id", "pct")
cell_BC_x = melt(as.matrix(bc_matrix))
cell_BC$UMI_count = as.vector(cell_BC_x$value)
cell_BC = cell_BC %>% group_by(cell_id) %>% arrange(desc(pct), .by_group = TRUE)
cell_BC$pct_order = rep(c(1:nrow(bc_matrix)), times = ncol(bc_matrix))

### pct >= 30 and UMI_count >= 10
cell_BC_sub = subset(cell_BC, pct >= 30 & UMI_count >= 10)
print(length(unique(cell_BC_sub$cell_id))/length(unique(cell_BC$cell_id)))
print(sum(table(cell_BC_sub$cell_id) == 1)/length(unique(cell_BC$cell_id)))
### 77.7% cells are retained, and 69.2% cells have unique BC

### calculating Peasrson correlation on the UMI matrix
BC_list_num = cell_BC_sub %>% group_by(BC) %>% tally() %>% arrange(desc(n)) %>% 
    mutate(log2_n = log2(n)) %>% as.data.frame()
BC_list_num$BC_id = paste0("BC_",1:nrow(BC_list_num))
BC_list = as.vector(BC_list_num$BC)

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

### identifying clones with MOI > 1
dat_x = dat %>% filter(A != B, corr > 0.1) %>% arrange(desc(corr)) 

### identifying distinct clonotypes

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


### merging clones with MOI > 1

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

saveRDS(cell_BC_assigned, paste0(work_path, "/analysis/monoclonal_EB_TF_screen_sci/cell_to_mEB/cell_clonotype_rep1.rds"))
print(nrow(cell_BC_assigned)/sum(pd$RT_group %in% c("GEx_mEB_ebTF_SR_filtered_handpool_1to10","GEx_mEB_ebTF_SR_gravity_handpool_1to10")))
### Finally, 17.3% of cells with transcriptome have been assigned to a distinct clonotype in rep1
### 49 clonotypes
### mean 88, median 65, sd 58, ranging 50 - 302 cells per clonotypes


############################################
### assigning gRNAs for individual clonotype

### gRNA_BC look up table
gRNA_BC = read.table("/net/shendure/vol8/projects/cxqiu/nobackup/work/sam_tf/processing/monoclonal_EB_TF_screen_sci/gRNA_BC_look_up_table.txt", as.is=T)
colnames(gRNA_BC) = c("gRNA", "BC")

reverse_complement_list <- function(dna_sequences) {
    complement_bases <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
    reverse_complement_single <- function(sequence) {
        split_sequence <- unlist(strsplit(sequence, ""))
        complement_sequence <- complement_bases[split_sequence]
        reverse_complement_sequence <- rev(complement_sequence)
        paste(reverse_complement_sequence, collapse = "")
    }
    lapply(dna_sequences, reverse_complement_single)
}
gRNA_BC$BC_rc = unlist(reverse_complement_list(gRNA_BC$BC))


cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen_sci/cell_to_mEB/cell_clonotype_rep1.rds"))
cell_BC_assigned = cell_BC_assigned[,c("BC", "clonotype_id", "cell_id")]
cell_BC_assigned = tidyr::separate_rows(cell_BC_assigned, BC, sep = ",")

library(stringr)
cell_BC_assigned$BC_rc = str_sub(as.vector(cell_BC_assigned$BC), -9, -2)
cell_BC_assigned = cell_BC_assigned %>% left_join(gRNA_BC[,c("gRNA","BC_rc")], by = "BC_rc")

cell_BC_assigned$TF = unlist(lapply(as.vector(cell_BC_assigned$gRNA), function(x) strsplit(x,"[.]")[[1]][1])) 

cell_BC_assigned = cell_BC_assigned %>% select(clonotype_id, BC, gRNA, TF, sample = cell_id) %>%
    left_join(pd %>% select(sample, celltype))


x = cell_BC_assigned %>% select(clonotype_id, TF) %>% unique() %>%
    group_by(TF) %>% tally()



#########################
### box plot

target_i = "Hand1"
celltype_i = "Cardiomyocytes"

df_x = cell_BC_assigned %>% group_by(TF, clonotype_id, celltype) %>% tally() %>%
    left_join(cell_BC_assigned %>% group_by(clonotype_id) %>% tally() %>% rename(total_n=n), by = "clonotype_id") %>%
    mutate(pct = 100*n/total_n)

try(ggplot(df_x %>% filter(TF %in% c(target_i, "NONTARGETING"), celltype == celltype_i), aes(TF, pct, fill = TF)) + geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2) + 
        labs(x="", y="% of cells", title = celltype_i) +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        scale_fill_brewer(palette="Paired") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen_sci/cell_to_mEB/rep1_boxplot_", target_i, "_", gsub(" ","_",celltype_i), ".pdf"),
               height  = 5, 
               width = 3), silent = T)

x_1 = as.vector(df_x$pct[df_x$TF == target_i & df_x$celltype == celltype_i])
x_2 = as.vector(df_x$pct[df_x$TF == "NONTARGETING" & df_x$celltype == celltype_i])
wilcox.test(x_1, x_2) ### p = 0.063




target_i = "Hand1"
celltype_i = "Lateral plate mesoderm"

df_x = cell_BC_assigned %>% group_by(TF, clonotype_id, celltype) %>% tally() %>%
    left_join(cell_BC_assigned %>% group_by(clonotype_id) %>% tally() %>% rename(total_n=n), by = "clonotype_id") %>%
    mutate(pct = 100*n/total_n)

try(ggplot(df_x %>% filter(TF %in% c(target_i, "NONTARGETING"), celltype == celltype_i), aes(TF, pct, fill = TF)) + geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2) + 
        labs(x="", y="% of cells", title = celltype_i) +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        scale_fill_brewer(palette="Paired") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
        ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen_sci/cell_to_mEB/rep1_boxplot_", target_i, "_", gsub(" ","_",celltype_i), ".pdf"),
               height  = 5, 
               width = 3), silent = T)


x_1 = as.vector(df_x$pct[df_x$TF == target_i & df_x$celltype == celltype_i])
x_2 = as.vector(df_x$pct[df_x$TF == "NONTARGETING" & df_x$celltype == celltype_i])
wilcox.test(x_1, x_2) ### p = 0.029


