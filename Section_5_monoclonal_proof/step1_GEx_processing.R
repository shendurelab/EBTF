
#######################################
### Processing monoclonal_proof dataset
### Chengxiang Qiu
### Feb-20, 2025

### processing monoclonal_proof dataset, including two technical replicates (rep1_1 and rep1_2)
### Of note, gRNA and transcriptome were sequenced together from 10x, while bc was sequenced separately

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

run_list = paste0("E2E_rep", c("1_1", "1_2"))

###########################
### Step-1: read 10x output

for(kk in c(1:2)){
    run_id = run_list[kk]
    print(run_id)
    
    count_matrix_all = Read10X(paste0(work_path, "/data/monoclonal_EB_proof/", run_id, "/outs/raw_feature_bc_matrix"), 
                               gene.column = 1,
                               strip.suffix = T)
    
    count_matrix = count_matrix_all[["Gene Expression"]]
    
    colnames(count_matrix) = paste0("monoclonal_EB_proof_", run_id, "_", colnames(count_matrix))
    
    count_matrix_nonzero = count_matrix
    count_matrix_nonzero@x[count_matrix_nonzero@x >= 1] = 1
    
    pd = data.frame(cell_id = colnames(count_matrix),
                    UMI_count = colSums(count_matrix),
                    gene_count = colSums(count_matrix_nonzero))
    rownames(pd) = as.vector(pd$cell_id)
    pd$experiment_id = paste0("monoclonal_EB_proof_", run_id)
    
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
    
    Matrix::writeMM(t(count_matrix_filter), paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/gene_count.mtx"))
    write.csv(pd_filter, paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/df_cell.csv"))
    
    saveRDS(pd_filter, paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/df_cell.rds"))
    
}

###########################
### Step-2: detect doublets

### running scrublet to identify doublets
### python detect_doublet.py

################################
### Step-3: plot summary results

for(run_id in run_list){
    print(run_id)
    
    pd_filter = readRDS(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/df_cell.rds"))
    scrublet_score = read.csv(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/doublet_scores_observed_cells.csv"), header=F)
    pd_filter$doublet_score = as.vector(scrublet_score$V1)
    saveRDS(pd_filter, paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/df_cell.rds"))
    
    simulated_scrublet_score = read.csv(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/doublet_scores_simulated_doublets.csv"), header=F)
    df = data.frame(score = c(as.vector(scrublet_score$V1), as.vector(simulated_scrublet_score$V1)),
                    group = rep(c("observed", "simulated"), times = c(nrow(scrublet_score), nrow(simulated_scrublet_score))))
    p = ggplot(df, aes(score, fill = group)) + geom_histogram(alpha=0.6, binwidth = 0.02)
    pdf(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/hist_doublet_score.pdf"))
    print(p)
    dev.off()
    
    pdf(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/hist_log2_umi_orig.pdf"))
    print(hist(log2(pd$UMI_count + 1), 100))
    dev.off()
    
    pdf(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/hist_log2_umi.pdf"))
    print(hist(pd_filter$log2_umi, 100))
    dev.off()
    
    pdf(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/hist_MT_pct.pdf"))
    print(hist(pd_filter$MT_pct, 100))
    dev.off()
    
    pdf(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/hist_RIBO_pct.pdf"))
    print(hist(pd_filter$RIBO_pct, 100))
    dev.off()
    
}


##########################
### Step-4:filtering cells 


for(run_id in run_list){    
    print(run_id)
        
    count_matrix_all = Read10X(paste0(work_path, "/data/monoclonal_EB_proof/", run_id, "/outs/raw_feature_bc_matrix"), 
                           gene.column = 1,
                           strip.suffix = T)
    count_matrix = count_matrix_all[["Gene Expression"]]
    
    colnames(count_matrix) = paste0("monoclonal_EB_proof_", run_id, "_", colnames(count_matrix))
    
    pd_filter = readRDS(paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/df_cell.rds"))
    count_matrix_filter = count_matrix[rownames(count_matrix) %in% as.vector(mouse_gene_sub$gene_ID),]
    
    umi_cutoff = 1000
    MT_pct_cutoff_up = 10
    RIBO_pct_cutoff = 40
    
    pd_filter = pd_filter[pd_filter$UMI_count >= umi_cutoff &
                          pd_filter$UMI_count <= quantile(pd_filter$UMI_count, 0.995) &
                          pd_filter$doublet_score <= 0.2 &
                          pd_filter$MT_pct <= MT_pct_cutoff_up &
                          pd_filter$RIBO_pct <= RIBO_pct_cutoff,]
    
    count_matrix_filter = count_matrix_filter[,as.vector(pd_filter$cell_id)]
    
    obj = CreateSeuratObject(count_matrix_filter, meta.data = pd_filter)
    saveRDS(obj, paste0(work_path, "/processing/monoclonal_EB_proof/", run_id, "/obj_", "monoclonal_EB_proof_", run_id, ".rds"))
    
}


run_list=(rep1_1 rep1_2 rep2_1 rep2_2)
for i in "${run_list[@]}"; do echo "$i"; done
for i in "${run_list[@]}";do rm ./E2E_"$i"/df_cell.csv ./E2E_"$i"/gene_count.mtx;done
for i in "${run_list[@]}";do mkdir -p ./E2E_"$i"/doublet_removing;done
for i in "${run_list[@]}";do mv ./E2E_"$i"/*.pdf ./E2E_"$i"/*.csv ./E2E_"$i"/doublet_removing;done







