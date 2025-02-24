
#######################################################
### Processing monoclonal_screen dataset (sci-RNA-seq3)
### Chengxiang Qiu
### Feb-20, 2025

####################################
### Step-1, creating the data matrix

source("help_script.R")

batch_num = 1

### how many reads in this experiment (after UMI attach)
read_num_fastq = NULL
for(cnt in 1:batch_num){
    print(cnt)
    read_num_cnt = read.table(paste0(data_path, "/nobackup/output_", cnt, "/read_num_UMI_attach.txt"),as.is=T)
    read_num_fastq = rbind(read_num_fastq, read_num_cnt)
}
print(sum(read_num_fastq$V1))
### 104,817,124

### summary the duplication rate
read_num = NULL
for(cnt in 1:batch_num){
  print(cnt)
  read_num_cnt = read.table(paste0(data_path, "/nobackup/output_", cnt, "/read_num.txt"),as.is=T)
  read_num = rbind(read_num, read_num_cnt)
}

print(summary(1 - read_num$V2/read_num$V1)) 
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
### 0.06642 0.19272 0.20152 0.19748 0.21094 0.23987

print(sum(read_num$V1))
### 68,786,065

print(sum(read_num$V2))
### 54,791,551

### saving data and performing doublets removing on individual batches if necessary

df_cell_merge = NULL

for(cnt in 1:batch_num){
  print(cnt)
    
  load(paste0(data_path, "/nobackup/output_", cnt, "/report/sci_summary.RData"))

  print(sum(rownames(gene_count) != df_gene$gene_id))
  print(sum(colnames(gene_count) != df_cell$sample))
  
  rownames(df_gene) = unlist(lapply(rownames(df_gene), function(x) strsplit(x,"[.]")[[1]][1]))
  rownames(gene_count) = unlist(lapply(rownames(gene_count), function(x) strsplit(x,"[.]")[[1]][1]))
  df_gene$gene_id = unlist(lapply(as.vector(df_gene$gene_id), function(x) strsplit(x,"[.]")[[1]][1]))
  
  print(dim(gene_count))
  
  ### read RT_barcode and bbi-sample-sheet to extract sample name (or RT_group) for each individual cell
  RT_plate = read.table(paste0(data_path, "/rt.txt"), as.is=T)
  colnames(RT_plate) = c("RT_barcode_well", "RT_barcode_sequence")
  
  LIG_plate = read.table(paste0(data_path, "/ligation.txt"), as.is=T)
  colnames(LIG_plate) = c("LIG_barcode_well", "LIG_barcode_sequence")
  
  PCR_plate = read.table(paste0(data_path, "/pair_PCR_idx_GEx_seq072_v3_20240705.txt"), header=T, as.is=T, sep=",")
  colnames(PCR_plate) = c("PCR_barcode_well", "PCR_p5_sequence", "PCR_p7_sequence")
  
  RT_barcode = read.table(paste0(data_path, "/bbi_RT_sample_sheet_GEx_seq072_20240618.txt"), as.is=T, header=F, sep=",")
  names(RT_barcode) = c("RT_barcode_well", "RT_group", "RT_species")
  RT_barcode = RT_barcode %>% left_join(RT_plate, by = "RT_barcode_well") %>% as.data.frame()
  
  df_cell$RT_barcode_sequence = str_sub(unlist(lapply(as.vector(df_cell$sample), function(x) strsplit(x,"[.]")[[1]][2])), -10, -1)
  df_cell$LIG_barcode_sequence = str_sub(unlist(lapply(as.vector(df_cell$sample), function(x) strsplit(x,"[.]")[[1]][2])), 1, -11)
  df_cell$PCR_barcode_well = unlist(lapply(as.vector(df_cell$sample), function(x) strsplit(x,"[.]")[[1]][1]))
  df_cell$PCR_barcode_well = sub("_[^_]*$", "", df_cell$PCR_barcode_well)
  
  df_cell = df_cell %>%
    left_join(RT_barcode, by = "RT_barcode_sequence") %>%
    left_join(LIG_plate, by = "LIG_barcode_sequence") %>%
    left_join(PCR_plate, by = "PCR_barcode_well")
  print(sum(is.na(df_cell$RT_group)))

  df_cell = df_cell[,c("sample", "unmatched_rate", "all_exon", "all_intron",
                       "RT_barcode_well", "RT_barcode_sequence", "RT_group",
                       "LIG_barcode_well", "LIG_barcode_sequence",
                       "PCR_barcode_well", "PCR_p5_sequence", "PCR_p7_sequence")] 
  
  sum(is.na(df_cell$RT_group))
  sum(is.na(df_cell$LIG_barcode_well))
  sum(is.na(df_cell$PCR_p5_sequence))
  
  df_cell$UMI_count = Matrix::colSums(gene_count)
  gene_count_copy = gene_count
  gene_count_copy@x[gene_count_copy@x > 0] = 1
  df_cell$gene_count = Matrix::colSums(gene_count_copy)
  
  mouse_gene = read.table("/net/gs/vol1/home/cxqiu/work/tome/code/mouse.v12.geneID.txt", header=T, sep="\t", as.is=T)
  rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
  mouse_gene = mouse_gene[as.vector(df_gene$gene_id),]
  df_gene$chr = as.vector(mouse_gene$chr)
  names(df_gene) = c("gene_id", "gene_type", "gene_short_name", "chr")
  gene_keep = df_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y"))
  
  keep = df_cell$UMI_count >= 100 & 
    df_cell$gene_count >= 100 & 
    df_cell$unmatched_rate < 0.4 &
    df_cell$RT_group != "GEx_mEB_ebMPRA_day22"
  df_cell = df_cell[keep,]
  df_gene = df_gene[gene_keep,]
  gene_count = gene_count[gene_keep, keep]
  rownames(gene_count) = as.vector(df_gene$gene_id)
  colnames(gene_count) = as.vector(df_cell$sample)
  rownames(df_cell) = as.vector(df_cell$sample)
  
  print(dim(gene_count))
  
  saveRDS(gene_count, paste0(work_path, "/gene_count_", cnt, ".rds"))

  df_cell_merge = rbind(df_cell_merge, df_cell)  
}

print(nrow(df_cell_merge))
### nrow(df_cell) = 54,356

print(median(df_cell_merge$UMI_count))
### 555

print(median(df_cell_merge$gene_count))
### 465

saveRDS(df_cell_merge, paste0(work_path, "/df_cell.rds"))


###############################
### Step-2: detecting scrublets

### If the dataset is small, we don't need to split it into multiple batches for running scrublet

df_cell_merge = readRDS(paste0(work_path, "/df_cell.rds"))
gene_count_merge = NULL
for(cnt in 1:batch_num){
    gene_count = readRDS(paste0(work_path, "/gene_count_", cnt, ".rds"))
    gene_count_merge = cbind(gene_count_merge, gene_count)
}
print(sum(colnames(gene_count_merge) == rownames(df_cell_merge)))

writeMM(t(gene_count_merge), paste0(work_path, "/gene_count_1.mtx"))
write.csv(df_cell_merge, paste0(work_path, "/df_cell_1.csv"))


### run scrublet using python to detect doublets
### python detect_scrublet.py

df_cell = readRDS(paste0(work_path, "/df_cell.rds"))
batch_num = 1

doublet_scores_observed_cells = NULL
doublet_scores_simulated_doublets = NULL
df = NULL
for(i in 1:batch_num){
    print(i)
    df_i = read.csv(paste0(work_path, "/df_cell_", i, ".csv"), header=T, row.names=1, as.is=T)
    
    doublet_scores_observed_cells_i = read.csv(paste0(work_path, "/doublet_scores_observed_cells.csv"), header=F)
    df_i$doublet_score = as.vector(doublet_scores_observed_cells_i$V1)
    df = rbind(df, df_i)
    
    doublet_scores_simulated_doublets_i = read.csv(paste0(work_path, "/doublet_scores_simulated_doublets.csv"), header=F)
    doublet_scores_simulated_doublets = rbind(doublet_scores_simulated_doublets, doublet_scores_simulated_doublets_i)
}

df = df[rownames(df_cell),]
print(sum(rownames(df) != rownames(df_cell)))

df_cell$doublet_score = as.vector(df$doublet_score)
df_cell$detected_doublets = df_cell$doublet_score > 0.2

### sum(df_cell$detected_doublets)/nrow(df_cell) = 0.03944367

###############################################################
### checking if sub-clusters include over 15% doublet cells ###
###############################################################

global = read.csv(paste0(work_path, "/doublet_cluster/global.csv"), header=T)
main_cluster_list = sort(as.vector(unique(global$louvain)))

res = NULL

for(i in 1:length(main_cluster_list)){
    print(paste0(i, "/", length(main_cluster_list)))
    dat = read.csv(paste0(work_path, "/doublet_cluster/adata.obs.louvain_", (i-1), ".csv"), header=T)
    print(nrow(dat))
    dat$louvain = as.vector(paste0("cluster_", dat$louvain))
    dat = dat %>%
        left_join(df_cell[,c("sample", "detected_doublets", "doublet_score")], by = "sample")
    
    tmp2 = dat %>%
        group_by(louvain) %>%
        tally() %>%
        dplyr::rename(n_sum = n)
    
    tmp1 = dat %>%
        filter(detected_doublets == "TRUE") %>%
        group_by(louvain) %>%
        tally() %>%
        left_join(tmp2, by = "louvain") %>%
        mutate(frac = n/n_sum) %>%
        filter(frac > 0.15)
    
    dat$doublet_cluster = dat$louvain %in% as.vector(tmp1$louvain) 
    
    p1 = ggplot(dat, aes(umap_1, umap_2, color = louvain)) + geom_point() + theme(legend.position="none") 
    p2 = ggplot(dat, aes(umap_1, umap_2, color = doublet_cluster)) + geom_point()
    p3 = ggplot(dat, aes(umap_1, umap_2, color = detected_doublets)) + geom_point()
    p4 = ggplot(dat, aes(umap_1, umap_2, color = doublet_score)) + geom_point() + scale_color_viridis(option = "plasma")
    p5 = ggplot(dat, aes(umap_1, umap_2, color = UMI_count)) + geom_point() + scale_color_viridis(option = "plasma")
    p6 = ggplot(dat, aes(umap_1, umap_2, color = gene_count)) + geom_point() + scale_color_viridis(option = "plasma")
    
    pdf(paste0(work_path, "/doublet_cluster/adata.obs.louvain_", (i-1), ".pdf"), 12, 18)
    grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3, ncol=2) 
    dev.off()
    
    dat[,!colnames(dat) %in% c("umap_1", "umap_2")]
    dat$main_louvain = (i-1)
    
    res = rbind(res, dat)
}

rownames(res) = as.vector(res$sample)
res = res[rownames(df_cell),]
df_cell$doublet_cluster = res$doublet_cluster

### sum(df_cell$detected_doublets | df_cell$doublet_cluster) = 4086
### sum(df_cell$detected_doublets | df_cell$doublet_cluster)/nrow(df_cell) = 0.07517109
saveRDS(df_cell, paste0(work_path, "/df_cell.rds"))

### mv doublet_scores_observed_cells.csv doublet_cluster/
### mv doublet_scores_simulated_doublets.csv doublet_cluster/
### rm df_cell_1.csv gene_count_1.mtx




#####################################################
### Step-3: checking more cutoffs for filtering cells

pd = readRDS(paste0(work_path, "/df_cell.rds"))
count = readRDS(paste0(work_path, "/gene_count_1.rds"))
fd = mouse_gene[rownames(count),]
fd$gene_id = fd$gene_ID

pd$log2_umi = log2(pd$UMI_count)
pd$EXON_pct = 100 * pd$all_exon / (pd$all_exon + pd$all_intron)
print(nrow(pd))
### n = 54,356 cells

### calculate MT_pct and Ribo_pct per cell
gene = fd
MT_gene = as.vector(gene[grep("^mt-",gene$gene_short_name),]$gene_id)
Rpl_gene = as.vector(gene[grep("^Rpl",gene$gene_short_name),]$gene_id)
Mrpl_gene = as.vector(gene[grep("^Mrpl",gene$gene_short_name),]$gene_id)
Rps_gene = as.vector(gene[grep("^Rps",gene$gene_short_name),]$gene_id)
Mrps_gene = as.vector(gene[grep("^Mrps",gene$gene_short_name),]$gene_id)
RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)

print(sum(colnames(count) != rownames(pd)))
pd$MT_pct = 100 * Matrix::colSums(count[gene$gene_id %in% MT_gene, ])/Matrix::colSums(count)
pd$RIBO_pct = 100 * Matrix::colSums(count[gene$gene_id %in% RIBO_gene, ])/Matrix::colSums(count)

print(sum(pd$detected_doublets | pd$doublet_cluster)/nrow(pd))
pd = pd[!(pd$detected_doublets | pd$doublet_cluster),]
pd$detected_doublets = pd$doublet_cluster = NULL
print(nrow(pd))
saveRDS(pd, paste0(work_path, "/pd.rds"))

pd_sub = pd %>% filter(MT_pct < 10, RIBO_pct < 10, doublet_score < 0.15)
print(nrow(pd_sub)/nrow(pd))

x1 = quantile(pd_sub$log2_umi, 0.025)
x2 = quantile(pd_sub$log2_umi, 0.975)
pd_sub = pd_sub %>% filter(log2_umi >= x1, log2_umi <= x2)
print(nrow(pd_sub)/nrow(pd))

pd = pd[pd$sample %in% pd_sub$sample,]

### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
fd = fd[(fd$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & fd$chr %in% paste0("chr", c(1:19, "M", "X", "Y")),]

print(sum(!colnames(count) %in% rownames(pd)))
count = count[rownames(count) %in% rownames(fd),rownames(pd)]

obj = CreateSeuratObject(count, meta.data = pd)
print(dim(obj))

saveRDS(obj, paste0(work_path, "/obj.rds"))


###############################
### check cutoff of UMI #######
###############################

pd = readRDS("pd.rds")
print(dim(pd))

p1 = ggplot(pd, aes(log2_umi)) + geom_histogram(binwidth = 0.1) 

p2 = ggplot(pd, aes(RIBO_pct)) + geom_histogram(binwidth = 0.1)

p3 = ggplot(pd, aes(MT_pct)) + geom_histogram(binwidth = 0.1) 

p4 = ggplot(pd, aes(EXON_pct)) + geom_histogram(binwidth = 0.1) + geom_vline(xintercept = 85) 

pdf("quality_summary.pdf",6,5)
grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2) 
dev.off()

pd_sub = pd %>% filter(MT_pct < 10, RIBO_pct < 10, doublet_score < 0.15)
print(nrow(pd_sub)/nrow(pd))

x1 = quantile(pd_sub$log2_umi, 0.025)
x2 = quantile(pd_sub$log2_umi, 0.975)

hist(pd_sub$log2_umi, 500); abline(v = x1); abline(v = x2)
print(sum(pd_sub$log2_umi >= x1 & pd_sub$log2_umi <= x2)/nrow(pd_sub))


##################################################
### Step-4: performing regular dimension reduction

obj = readRDS(paste0(work_path, "/obj.rds"))

count = GetAssayData(obj, slot = "counts")
df_gene_sub = mouse_gene[rownames(count),]
df_gene_sub = df_gene_sub[!df_gene_sub$chr %in% c("chrX", "chrY"),]
count = count[rownames(count) %in% rownames(df_gene_sub),]
obj = CreateSeuratObject(count, meta.data = data.frame(obj[[]]))

obj$group = "monoclonal_EB_TF_screen_sci"
obj_processed = doClusterSeurat(obj)
obj_processed = FindClusters(object = obj_processed, resolution = 2)
obj_processed = FindClusters(object = obj_processed, resolution = 5)

obj_processed$UMAP_2d_1 = Embeddings(obj_processed, reduction = "umap")[,1]
obj_processed$UMAP_2d_2 = Embeddings(obj_processed, reduction = "umap")[,2]
obj_processed = RunUMAP(object = obj_processed, 
                        reduction = "pca", 
                        dims = 1:30, 
                        min.dist = 0.3, 
                        n.components = 3)

obj_processed$UMAP_1 = Embeddings(obj_processed, reduction = "umap")[,1]
obj_processed$UMAP_2 = Embeddings(obj_processed, reduction = "umap")[,2]
obj_processed$UMAP_3 = Embeddings(obj_processed, reduction = "umap")[,3]

saveRDS(obj_processed, paste0(work_path, "/obj_processed.rds"))
saveRDS(data.frame(obj_processed[[]]), paste0(work_path, "/obj_processed_pd.rds"))


######################################
### creating a cds object for UMAP vis
genes_include = VariableFeatures(obj_processed)
cds = doObjectTransform(obj_processed, transform_to = "monocle3")
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_2d_1", "UMAP_2d_2")])
saveRDS(cds, paste0(work_path, "/cds.rds"))


################################
### Step-5: cell type annotation

pd = readRDS(paste0(work_path, "/obj_processed_pd.rds"))

anno = rep(NA, nrow(pd))
anno[pd$RNA_snn_res.1 %in% c(8)] = "Primordial germ cells"
anno[pd$RNA_snn_res.1 %in% c(13)] = "Epiblast"
anno[pd$RNA_snn_res.1 %in% c(11,18)] = "Primitive streak"
anno[pd$RNA_snn_res.1 %in% c(16)] = "Neuroectoderm"
anno[pd$RNA_snn_res.1 %in% c(12)] = "Surface ectoderm"
anno[pd$RNA_snn_res.1 %in% c(5,15)] = "Endoderm"
anno[pd$RNA_snn_res.1 %in% c(0,1,2)] = "Lateral plate mesoderm"
anno[pd$RNA_snn_res.1 %in% c(3,9)] = "Cardiomyocytes"
anno[pd$RNA_snn_res.1 %in% c(4)] = "Paraxial mesoderm" ### Meox1, Meox2, Tbx1, Pax3
anno[pd$RNA_snn_res.1 %in% c(7,10)] = "Endothelial cells"
anno[pd$RNA_snn_res.1 %in% c(6)] = "Erythroid cells"
anno[pd$RNA_snn_res.1 %in% c(14,17)] = "Blood progenitors"

anno[pd$RNA_snn_res.2 %in% c(5,10,13,14,19,21)] = "Cardiomyocytes"

pd$celltype = as.vector(anno)
saveRDS(pd, paste0(work_path, "/obj_processed_pd.rds"))


try(ggplot(pd) +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.12, color = "black") +
        geom_point(aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.08) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/monoclonal_EB_TF_screen_sci_celltype.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

pd_sub = pd %>% group_by(celltype) %>% summarize(UMAP_1_mean = mean(UMAP_2d_1), UMAP_2_mean = mean(UMAP_2d_2))
try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.12, color = "black") +
        geom_point(data = pd,
                   aes(x = UMAP_2d_1, y = UMAP_2d_2, color = celltype), size=0.08) +
        ggrepel::geom_text_repel(data = pd_sub, aes(x = UMAP_1_mean, y = UMAP_2_mean, label = celltype), color = "black", size = 2.5, family = "Arial") +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/plot/monoclonal_EB_TF_screen_sci_celltype_label.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


