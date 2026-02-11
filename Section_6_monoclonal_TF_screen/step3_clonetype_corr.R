
### Reviewer-1: last question

### I am confused by what the UMAPs in Fig 4c and 4d are showing. 
### My understanding from the figure legend and text is that these show 
### dimensionality reduction of the matrix of cell barcodes x clonotype barcode 
### UMIs (rather than of cell barcodes x transcriptomic UMIs). In principle, 
### each cell should have only 1, and at most several integrated barcodes. 
### A UMAP seems unnecessary to represent this data (and interpreting 85% 
### 'distance' on a UMAP is shaky, given that UMAP distances famously have no 
### physical interpretation). Is there a more straightforward way of validating 
### the correspondence of clonotype barcodes and perturbation conditions? For 
### example, the pairwise correlations between different clonotype barcodes?



######################
### Rep-1 (arrayed)

source("~/work/scripts/utils.R")
library(igraph)
mouse_gene <- read.table("~/work/tome/code/mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf"

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))

bc_matrix = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/BC_matrix_rep1.rds"))
gRNA_matrix = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/gRNA_matrix_rep1.rds"))

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
### 85.6% cells are retained, and 73.9% cells have unique BC

### calculating Peasrson correlation on the UMI matrix

BC_list_num = cell_BC_sub %>% group_by(BC) %>% tally() %>% arrange(desc(n)) %>% 
    mutate(log2_n = log2(n)) %>% as.data.frame()

BC_list_num = BC_list_num[1:50,]

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
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/BC_cell_num_rep1_top50.pdf"),
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
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/heatmap_UMI_correlation_BC_rep1_top50.pdf"),
           height  = 3, 
           width = 4)







######################
### Rep-2 (pooled)

source("~/work/scripts/utils.R")
library(igraph)
mouse_gene <- read.table("~/work/tome/code/mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf"

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))

bc_matrix = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/BC_matrix_rep2.rds"))
gRNA_matrix = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/gRNA_matrix_rep2.rds"))

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
### 74.8% cells are retained, and 69.3% cells have unique BC

### calculating Peasrson correlation on the UMI matrix

BC_list_num = cell_BC_sub %>% group_by(BC) %>% tally() %>% arrange(desc(n)) %>% 
    mutate(log2_n = log2(n)) %>% as.data.frame()

BC_list_num = BC_list_num[1:50,]

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
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/BC_cell_num_rep2_top50.pdf"),
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
    ggsave(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/heatmap_UMI_correlation_BC_rep2_top50.pdf"),
           height  = 3, 
           width = 4)


