
#############################################
### Identifying cell-type composition changes
### Chengxiang Qiu
### cxqiu@uw.edu

### Please read the tutorial of Hooke package:
### https://cole-trapnell-lab.github.io/projects/hooke/
### https://cole-trapnell-lab.github.io/papers/duran-genetic-inferance/

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

library(PLNmodels, lib.loc = "/net/gs/vol1/home/cxqiu/R/x86_64-pc-linux-gnu-library/4.3")
library(hooke)

#################################################################
### Step-1: Assigning gRNAs and targets for individual colonotype
### (taking the "arrayed" monoclonal EBs dataset as an example)

cell_gRNA = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep1.rds"))
cell_gRNA$target = unlist(lapply(as.vector(cell_gRNA$gRNA), function(x) strsplit(x,"[_]")[[1]][1])) 

cell_gRNA_uniq = cell_gRNA %>% filter(pct >= 70) %>% arrange(clonotype_id)
cell_gRNA_mult = cell_gRNA %>% filter(!clonotype_id %in% as.vector(cell_gRNA_uniq$clonotype_id), order %in% c(1,2))

### out of 202 clonotypes, 187 are assigned with uique gRNAs, 15 are assigned with double gRNAs
cell_gRNA = rbind(cell_gRNA_uniq, cell_gRNA_mult)
cell_target = cell_gRNA %>% select(clonotype_id, target) %>% unique()

### we further set a new cutoff of cell number for identifying EBs (cell # = 100)
cell_BC_assigned = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_rep1.rds"))
cell_BC_assigned_num = cell_BC_assigned %>% group_by(clonotype_id) %>% tally() %>% filter(n >= 100)
cell_target = cell_target[cell_target$clonotype_id %in% as.vector(cell_BC_assigned_num$clonotype_id),]
cell_BC_assigned = cell_BC_assigned[cell_BC_assigned$clonotype_id %in% as.vector(cell_BC_assigned_num$clonotype_id),]


#####################################################################
### Step-2: performing hooker analysis to identify changed cell types

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
obj = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))

count = GetAssayData(obj, slot = "counts")[,as.vector(cell_BC_assigned$cell_id)]
pd_sub = pd[as.vector(cell_BC_assigned$cell_id),]
pd_sub$clonotype_id = as.vector(cell_BC_assigned$clonotype_id)
pd_sub = pd_sub %>% select(cell_id, experiment_id, UMAP_1, UMAP_2, celltype, clonotype_id) %>% 
    left_join(cell_target, by = "clonotype_id") %>% as.data.frame()
rownames(pd_sub) = as.vector(pd_sub$cell_id)

cds = new_cell_data_set(count,
                        cell_metadata = pd_sub)
cds = preprocess_cds(cds, use_genes = VariableFeatures(obj))
reducedDims(cds)$UMAP = as.matrix(pData(cds)[,c("UMAP_1", "UMAP_2")])

ccs = new_cell_count_set(cds, 
                         sample_group = "clonotype_id", 
                         cell_group = "celltype")

### Carm1        Gata6        Hand1         Hes5         Hhex         Lhx1
### 16           14            8           14           17           16
### Lmo2 NONTARGETING         Six3            T
### 11           19           18           21

ccm  = new_cell_count_model(ccs,
                            main_model_formula_str = "~target")

res = NULL
for (target_i in c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1")){
    print(target_i)
    cond_exp = estimate_abundances(ccm, tibble::tibble(target = target_i))
    cond_not_exp = estimate_abundances(ccm, tibble::tibble(target = "NONTARGETING"))
    
    cond_ne_v_e_tbl = compare_abundances(ccm, cond_not_exp, cond_exp)

    res = rbind(res, cond_ne_v_e_tbl %>% select(cell_group, target_x, target_y,
                                                delta_log_abund, delta_log_abund_se, delta_p_value, delta_q_value))
}
res$fdr = p.adjust(res$delta_p_value, method = "fdr")
res %>% filter(fdr < 0.1) %>% arrange(target_y, delta_log_abund) %>% as.data.frame()
saveRDS(res, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/hooke_result.rds"))

write.table(res, paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/hooke_result.txt"), row.names=F, quote=F, sep="\t")




