
#######################################################
### Processing mouse EBs in the pilot screen experiment
### Chengxiang Qiu
### Feb-20, 2025


###############################
### Step-1: processing 10x data

### rep1 ###
### 1 CRISPRcut unsorted; 1 CRISPRi unsorted
### 2 CRISPRcut sorted; 2 CRISPRi sorted

module load cellranger/7.2.0

run_id=EB21_CRISPRcut_unsorted
run_id=EB21_CRISPRi_unsorted

cellranger count \
--fastqs /net/shendure/vol8/projects/ajh24/ajh24/proj/2018eb_factor_screen/data/reads/2018_10_22_eb_pilot_screen/HWGC5BGX5/outs/fastq_path \
--nosecondary \
--localcores 8 \
--include-introns true \
--sample "$run_id" \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_GEx/rep1_"$run_id" \
--id "$run_id" \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0

run_id=EB21_CRISPRcut_sorted
run_id=EB21_CRISPRi_sorted

cellranger count \
--fastqs /net/shendure/vol8/projects/ajh24/ajh24/proj/2018eb_factor_screen/data/reads/2018_10_22_eb_pilot_screen/HWGC5BGX5/outs/fastq_path \
--nosecondary \
--localcores 8 \
--include-introns true \
--sample EB21_CRISPRcut_sorted_1,EB21_CRISPRcut_sorted_2 \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_GEx/rep1_"$run_id" \
--id "$run_id" \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0


### rep2 ###
### 3 CRISPRcut sorted; 3 CRISPRi sorted
### 2 Nextseq runs for deeper coverage

module load cellranger/7.2.0

cellranger count \
--fastqs /net/shendure/vol10/projects/Samuel/nobackup/10X/2019_EB21_pilot_screen/HYJYBGXB_cells_1/outs/fastq_path,/net/shendure/vol10/projects/Samuel/nobackup/10X/2019_EB21_pilot_screen/HJCFKBGXB_cells_2/outs/fastq_path \
--nosecondary \
--localcores 8 \
--include-introns true \
--sample EB21_CRISPRcut_sorted_1,EB21_CRISPRcut_sorted_2,EB21_CRISPRcut_sorted_3 \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_GEx/rep2_EB21_CRISPRcut_sorted \
--id EB21_CRISPRcut_sorted \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0



###############################################
### Step-2: processing pilot_screen_GEx dataset

### detecting doublets, filtering cells by UMI count, mito%, and ribo% if necessary

source("help_script.R")
mouse_gene <- read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

run_id = "rep2_EB21_CRISPRi_sorted"

count_matrix = Read10X(paste0(work_path, "/data/pilot_screen_GEx/", run_id, "/outs/raw_feature_bc_matrix"), 
                       gene.column = 1,
                       strip.suffix = T)
colnames(count_matrix) = paste0("pilot_screen_GEx_", run_id, "_", colnames(count_matrix))

count_matrix_nonzero = count_matrix
count_matrix_nonzero@x[count_matrix_nonzero@x >= 1] = 1

pd = data.frame(cell_id = colnames(count_matrix),
                UMI_count = colSums(count_matrix),
                gene_count = colSums(count_matrix_nonzero))
rownames(pd) = as.vector(pd$cell_id)
pd$experiment_id = paste0("pilot_screen_GEx_", run_id)

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]
#write.csv(mouse_gene_sub, paste0(work_path, "/processing/df_gene.csv"))

count_matrix_filter = count_matrix[rownames(count_matrix) %in% as.vector(mouse_gene_sub$gene_ID),
                                   pd$UMI_count >= 500 & pd$gene_count >= 250]
pd_filter = pd[pd$UMI_count >= 500 & pd$gene_count >= 250,]
pd_filter$log2_umi = log2(pd_filter$UMI_count)

MT_gene = as.vector(mouse_gene[grep("^mt-",mouse_gene$gene_short_name),]$gene_ID)
Rpl_gene = as.vector(mouse_gene[grep("^Rpl",mouse_gene$gene_short_name),]$gene_ID)
Mrpl_gene = as.vector(mouse_gene[grep("^Mrpl",mouse_gene$gene_short_name),]$gene_ID)
Rps_gene = as.vector(mouse_gene[grep("^Rps",mouse_gene$gene_short_name),]$gene_ID)
Mrps_gene = as.vector(mouse_gene[grep("^Mrps",mouse_gene$gene_short_name),]$gene_ID)
RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)

pd_filter$MT_pct = 100 * Matrix::colSums(count_matrix_filter[rownames(count_matrix_filter) %in% MT_gene,])/Matrix::colSums(count_matrix_filter)
pd_filter$RIBO_pct = 100 * Matrix::colSums(count_matrix_filter[rownames(count_matrix_filter) %in% RIBO_gene,])/Matrix::colSums(count_matrix_filter)

Matrix::writeMM(t(count_matrix_filter), paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/gene_count.matrix"))
write.csv(pd_filter, paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/df_cell.csv"))

### mkdir doublet_removing

### running scrublet to identify doublets
### python detect_scrublet.py
### rm df_cell.csv gene_count.matrix

scrublet_score = read.csv(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/doublet_removing/doublet_scores_observed_cells.csv"), header=F)
pd_filter$doublet_score = as.vector(scrublet_score$V1)
saveRDS(pd_filter, paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/doublet_removing/df_cell.rds"))

simulated_scrublet_score = read.csv(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/doublet_removing/doublet_scores_simulated_doublets.csv"), header=F)
df = data.frame(score = c(as.vector(scrublet_score$V1), as.vector(simulated_scrublet_score$V1)),
                group = rep(c("observed", "simulated"), times = c(nrow(scrublet_score), nrow(simulated_scrublet_score))))
p = ggplot(df, aes(score, fill = group)) + geom_histogram(alpha=0.6, binwidth = 0.02)
pdf(paste0(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/doublet_removing/doublet_score_hist.pdf")))
print(p)
dev.off()

### filtering cells 

if(run_id %in% c("rep1_EB21_CRISPRcut_unsorted", "rep1_EB21_CRISPRi_unsorted")){
    MT_cutoff = 10
} else {
    MT_cutoff = 20
}

pd_filter = pd_filter[pd_filter$UMI_count >= 1000 &
                      pd_filter$UMI_count <= quantile(pd_filter$UMI_count, 0.995) &
                      pd_filter$doublet_score  <= 0.2 &
                      pd_filter$MT_pct <= MT_cutoff,]

count_matrix_filter = count_matrix_filter[,as.vector(pd_filter$cell_id)]

obj = CreateSeuratObject(count_matrix_filter, meta.data = pd_filter)

saveRDS(obj, paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/obj_", "pilot_screen_GEx_", run_id, ".rds"))


