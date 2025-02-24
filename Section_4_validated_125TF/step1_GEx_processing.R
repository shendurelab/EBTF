
##########################################################
### Processing mouse EBs in the validated 125TF experiment
### Chengxiang Qiu
### Feb-20, 2025

##############################################
### Step-1: Processing data using 10x pipeline

#!/bin/bash
INPUTFILES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32)
INPUTFILENAME="${INPUTFILES[$SGE_TASK_ID - 1]}"
echo $INPUTFILENAME

module load cellranger/7.2.0

cellranger count \
--fastqs /net/shendure/vol9/nobackup/sam_tf/SD_SR_125TF_EB_10X/10X_test_done/outs/fastq_path \
--nosecondary \
--localcores 8 \
--include-introns true \
--sample EB_"$INPUTFILENAME" \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/EB21_CRISPRcut_125TF_GEx/EB_"$INPUTFILENAME" \
--id EB_"$INPUTFILENAME" \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0

### processing EB21_CRISPRcut_125TF_GEx dataset, including 30 lanes, each of which was pooled with multiple EBs
### detecting doublets, filtering cells by UMI count, mito%, and ribo% if necessary

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

run_list = paste0("EB_", c(1:32))

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])
run_id = run_list[kk]
print(run_id)

########################
### Step-2: reading data

count_matrix = Read10X(paste0(work_path, "/data/EB21_CRISPRcut_125TF_GEx/", run_id, "/outs/raw_feature_bc_matrix"), 
                       gene.column = 1,
                       strip.suffix = T)
colnames(count_matrix) = paste0("EB21_CRISPRcut_125TF_GEx_", run_id, "_", colnames(count_matrix))

count_matrix_nonzero = count_matrix
count_matrix_nonzero@x[count_matrix_nonzero@x >= 1] = 1

pd = data.frame(cell_id = colnames(count_matrix),
                UMI_count = colSums(count_matrix),
                gene_count = colSums(count_matrix_nonzero))
rownames(pd) = as.vector(pd$cell_id)
pd$experiment_id = paste0("EB21_CRISPRcut_125TF_GEx_", run_id)

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

Matrix::writeMM(t(count_matrix_filter), paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/gene_count.mtx"))
write.csv(pd_filter, paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/df_cell.csv"))

saveRDS(pd_filter, paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/df_cell.rds"))

##############################
### Step-3: detecting doublets

### running scrublet to identify doublets
### python detect_scrublet.py

########################################
### Step-4: Making quality summary plots

for(run_id in run_list){
    print(run_id)
    
    pd_filter = readRDS(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/df_cell.rds"))
    scrublet_score = read.csv(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/doublet_scores_observed_cells.csv"), header=F)
    pd_filter$doublet_score = as.vector(scrublet_score$V1)
    saveRDS(pd_filter, paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/df_cell.rds"))
    
    simulated_scrublet_score = read.csv(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/doublet_scores_simulated_doublets.csv"), header=F)
    df = data.frame(score = c(as.vector(scrublet_score$V1), as.vector(simulated_scrublet_score$V1)),
                    group = rep(c("observed", "simulated"), times = c(nrow(scrublet_score), nrow(simulated_scrublet_score))))
    p = ggplot(df, aes(score, fill = group)) + geom_histogram(alpha=0.6, binwidth = 0.02)
    pdf(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/hist_doublet_score.pdf"))
    print(p)
    dev.off()
    
    pdf(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/hist_log2_umi_orig.pdf"))
    print(hist(log2(pd$UMI_count + 1), 100))
    dev.off()
    
    pdf(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/hist_log2_umi.pdf"))
    print(hist(pd_filter$log2_umi, 100))
    dev.off()
    
    pdf(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/hist_MT_pct.pdf"))
    print(hist(pd_filter$MT_pct, 100))
    dev.off()
    
}


###########################
### Step-5: filtering cells


count_matrix = Read10X(paste0(work_path, "/data/EB21_CRISPRcut_125TF_GEx/", run_id, "/outs/raw_feature_bc_matrix"), 
                       gene.column = 1,
                       strip.suffix = T)
colnames(count_matrix) = paste0("EB21_CRISPRcut_125TF_GEx_", run_id, "_", colnames(count_matrix))

pd_filter = readRDS(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/df_cell.rds"))
count_matrix_filter = count_matrix[rownames(count_matrix) %in% as.vector(mouse_gene_sub$gene_ID),]

umi_cutoff = 2000
MT_pct_cutoff_up = 20
RIBO_pct_cutoff = 40

pd_filter = pd_filter[pd_filter$UMI_count >= umi_cutoff &
                      pd_filter$UMI_count <= quantile(pd_filter$UMI_count, 0.995) &
                      pd_filter$doublet_score <= 0.2 &
                      pd_filter$MT_pct <= MT_pct_cutoff_up &
                      pd_filter$RIBO_pct <= RIBO_pct_cutoff,]

count_matrix_filter = count_matrix_filter[,as.vector(pd_filter$cell_id)]

obj = CreateSeuratObject(count_matrix_filter, meta.data = pd_filter)

saveRDS(obj, paste0(work_path, "/processing/EB21_CRISPRcut_125TF_GEx/", run_id, "/obj_", "EB21_CRISPRcut_125TF_GEx_", run_id, ".rds"))




for i in $(seq 1 32);do rm ./EB_"$i"/df_cell.csv ./EB_"$i"/gene_count.mtx;done
for i in $(seq 1 32);do mkdir -p ./EB_"$i"/doublet_removing;done
for i in $(seq 1 32);do mv ./EB_"$i"/*.pdf ./EB_"$i"/*.csv ./EB_"$i"/doublet_removing;done



