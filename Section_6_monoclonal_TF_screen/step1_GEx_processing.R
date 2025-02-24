
########################################
### Processing monoclonal_screen dataset
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

### Of note, gRNA and transcriptome were sequenced together; 
### bc was sequenced separately, which has also been processed separately from here:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tf
### E2E_rep1p1_get_bc_v3_no_G_cleaned_UMI_20221213.txt

run_list = paste0("clonal_EBTF_", c(1:16))

####################################
### Step-1: read 10x pipeline output

for(kk in c(1:length(run_list))){
    run_id = run_list[kk]
    print(run_id)
    
    count_matrix_all = Read10X(paste0(work_path, "/data/monoclonal_EB_TF_screen/", run_id, "/outs/raw_feature_bc_matrix"), 
                               gene.column = 1,
                               strip.suffix = T)
    
    count_matrix = count_matrix_all[["Gene Expression"]]
    
    colnames(count_matrix) = paste0("monoclonal_EB_TF_screen_", run_id, "_", colnames(count_matrix))
    
    count_matrix_nonzero = count_matrix
    count_matrix_nonzero@x[count_matrix_nonzero@x >= 1] = 1
    
    pd = data.frame(cell_id = colnames(count_matrix),
                    UMI_count = colSums(count_matrix),
                    gene_count = colSums(count_matrix_nonzero))
    rownames(pd) = as.vector(pd$cell_id)
    pd$experiment_id = paste0("monoclonal_EB_TF_screen_", run_id)
    
    count_matrix_filter = count_matrix[rownames(count_matrix) %in% as.vector(mouse_gene_sub$gene_ID),
                                       pd$UMI_count >= 200 & pd$gene_count >= 100]
    pd_filter = pd[pd$UMI_count >= 200 & pd$gene_count >= 100,]
    pd_filter$log2_umi = log2(pd_filter$UMI_count)
    
    MT_gene = as.vector(mouse_gene[grep("^mt-",mouse_gene$gene_short_name),]$gene_ID)
    Rpl_gene = as.vector(mouse_gene[grep("^Rpl",mouse_gene$gene_short_name),]$gene_ID)
    Mrpl_gene = as.vector(mouse_gene[grep("^Mrpl",mouse_gene$gene_short_name),]$gene_ID)
    Rps_gene = as.vector(mouse_gene[grep("^Rps",mouse_gene$gene_short_name),]$gene_ID)
    Mrps_gene = as.vector(mouse_gene[grep("^Mrps",mouse_gene$gene_short_name),]$gene_ID)
    RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)
    
    pd_filter$MT_pct = 100 * Matrix::colSums(count_matrix_filter[rownames(count_matrix_filter) %in% MT_gene,])/Matrix::colSums(count_matrix_filter)
    pd_filter$RIBO_pct = 100 * Matrix::colSums(count_matrix_filter[rownames(count_matrix_filter) %in% RIBO_gene,])/Matrix::colSums(count_matrix_filter)
    
    Matrix::writeMM(t(count_matrix_filter), paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/gene_count.mtx"))
    write.csv(pd_filter, paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/df_cell.csv"))
    
    saveRDS(pd_filter, paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/df_cell.rds"))
    
}

###########################
### Step-2: detect doublets

### running scrublet to identify doublets
### python detect_doublet.py

####################################
### Step-3: checking quality summary

for(run_id in run_list){
    print(run_id)
    
    pd_filter = readRDS(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/df_cell.rds"))
    scrublet_score = read.csv(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/doublet_scores_observed_cells.csv"), header=F)
    pd_filter$doublet_score = as.vector(scrublet_score$V1)
    saveRDS(pd_filter, paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/df_cell.rds"))
    
    simulated_scrublet_score = read.csv(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/doublet_scores_simulated_doublets.csv"), header=F)
    df = data.frame(score = c(as.vector(scrublet_score$V1), as.vector(simulated_scrublet_score$V1)),
                    group = rep(c("observed", "simulated"), times = c(nrow(scrublet_score), nrow(simulated_scrublet_score))))
    p = ggplot(df, aes(score, fill = group)) + geom_histogram(alpha=0.6, binwidth = 0.02)
    pdf(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/hist_doublet_score.pdf"))
    print(p)
    dev.off()
    
    pdf(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/hist_log2_umi_orig.pdf"))
    print(hist(log2(pd$UMI_count + 1), 100))
    dev.off()
    
    pdf(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/hist_log2_umi.pdf"))
    print(hist(pd_filter$log2_umi, 100))
    dev.off()
    
    pdf(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/hist_MT_pct.pdf"))
    print(hist(pd_filter$MT_pct, 100))
    dev.off()
    
    pdf(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/hist_RIBO_pct.pdf"))
    print(hist(pd_filter$RIBO_pct, 100))
    dev.off()
    
}

###########################
### Step-4: filtering cells
 
for(run_id in run_list){    
    print(run_id)
        
    count_matrix_all = Read10X(paste0(work_path, "/data/monoclonal_EB_TF_screen/", run_id, "/outs/raw_feature_bc_matrix"), 
                           gene.column = 1,
                           strip.suffix = T)
    count_matrix = count_matrix_all[["Gene Expression"]]
    
    colnames(count_matrix) = paste0("monoclonal_EB_TF_screen_", run_id, "_", colnames(count_matrix))
    
    pd_filter = readRDS(paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/df_cell.rds"))
    count_matrix_filter = count_matrix[rownames(count_matrix) %in% as.vector(mouse_gene_sub$gene_ID),]
    
    if(run_id %in% paste0("clonal_EBTF_", c(1,2,3,4,5,8,9,10,12))){
        umi_cutoff = 1000
    } else if (run_id %in% paste0("clonal_EBTF_", c(6,7,11,13,15))){
        umi_cutoff = 750
    } else {
        umi_cutoff = 500
    }
    
    MT_pct_cutoff_up = 10
    RIBO_pct_cutoff = 50
    
    pd_filter = pd_filter[pd_filter$UMI_count >= umi_cutoff &
                          pd_filter$UMI_count <= quantile(pd_filter$UMI_count, 0.995) &
                          pd_filter$doublet_score <= 0.2 &
                          pd_filter$MT_pct <= MT_pct_cutoff_up &
                          pd_filter$RIBO_pct <= RIBO_pct_cutoff,]
    
    count_matrix_filter = count_matrix_filter[,as.vector(pd_filter$cell_id)]
    
    obj = CreateSeuratObject(count_matrix_filter, meta.data = pd_filter)
    saveRDS(obj, paste0(work_path, "/processing/monoclonal_EB_TF_screen/", run_id, "/obj_", "monoclonal_EB_TF_screen_", run_id, ".rds"))
    
}


run_list=(clonal_EBTF_1 clonal_EBTF_2 clonal_EBTF_3 clonal_EBTF_4 clonal_EBTF_5 clonal_EBTF_6 clonal_EBTF_7 clonal_EBTF_8 clonal_EBTF_9 clonal_EBTF_10 clonal_EBTF_11 clonal_EBTF_12 clonal_EBTF_13 clonal_EBTF_14 clonal_EBTF_15 clonal_EBTF_16)
for i in "${run_list[@]}";do echo "$i"; done
for i in "${run_list[@]}";do rm ./"$i"/df_cell.csv ./"$i"/gene_count.mtx;done
for i in "${run_list[@]}";do mkdir -p ./"$i"/doublet_removing;done
for i in "${run_list[@]}";do mv ./"$i"/*.pdf ./"$i"/*.csv ./"$i"/doublet_removing;done







