
######################################
### Power analysis based on simulation

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

library(tidyr)
library(forcats)
library(VGAM)
work_path = "Your_work_path"

################################################################################
### step-1, selecting subset of WT samples, used as reference for simulation ###
################################################################################

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
### 102,120 cells, including 12 cell types

### "Arrayed" condition assignments (6493 cells, 19 NTC samples)
gRNA_assign_1 = readRDS(paste0(work_path, 
                               "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep1_used_for_hooke.rds"))
gRNA_assign_1 = gRNA_assign_1[[1]] %>% left_join(gRNA_assign_1[[2]], by = "clonotype_id") %>%
    filter(target == "NONTARGETING", !is.na(target)) %>% mutate(clonotype = paste0("arrayed_", clonotype_id))

### "Pooled" condition assignments (1337 cells, 5 NTC samples)
gRNA_assign_2 = readRDS(paste0(work_path, 
                               "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep2_used_for_hooke.rds"))
gRNA_assign_2 = gRNA_assign_2[[1]] %>% left_join(gRNA_assign_2[[2]], by = "clonotype_id", relationship = "many-to-many") %>%
    filter(target == "NONTARGETING", !is.na(target)) %>% mutate(clonotype = paste0("pooled_", clonotype_id))

gRNA_assign = rbind(gRNA_assign_1[,c("cell_id", "clonotype")], 
                    gRNA_assign_2[,c("cell_id", "clonotype")])

### two NTC are larger than mean + 2*sd, so excluding them
### 4690 cells, 22 NTC samples
clonotype_cell_num = gRNA_assign %>% group_by(clonotype) %>% tally()
clonotype_cell_num = clonotype_cell_num %>%
    filter(n < mean(clonotype_cell_num$n) + 2*sd(clonotype_cell_num$n))
gRNA_assign = gRNA_assign[gRNA_assign$clonotype %in% clonotype_cell_num$clonotype,]

coldata_df = pd %>% filter(cell_id %in% gRNA_assign$cell_id) %>%
    select(cell_id, celltype) %>% left_join(gRNA_assign, by = "cell_id")

types = coldata_df %>% group_by(celltype) %>% tally() %>% 
    arrange(desc(n)) %>% pull(celltype)

prob = coldata_df %>% group_by(celltype) %>% tally() %>%
    mutate(percent = (n/sum(n))) %>% arrange(-percent) %>%
    filter(celltype %in% types) %>% pull(percent)

prob_save = coldata_df %>% group_by(celltype) %>% tally() %>%
    mutate(percent = (n/sum(n))) %>% arrange(-percent) %>%
    filter(celltype %in% types)

prob_x = coldata_df %>% group_by(celltype, clonotype) %>% tally() %>% dcast(celltype~clonotype, fill = 0)
rownames(prob_x) = prob_x[,1]
prob_x = as.matrix(prob_x[,-1])
y = apply(prob_x, 1, sd)/apply(prob_x, 1, mean)


###################################
### using Hooke analysis
### creating the cds and ccs object

library(PLNmodels, lib.loc = "/net/gs/vol1/home/cxqiu/R/x86_64-pc-linux-gnu-library/4.3")
library(hooke)

obj = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))
pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
count = GetAssayData(obj, slot = "counts")
keep = pd$cell_id %in% gRNA_assign$cell_id
pd_sub = pd[keep,]
count_sub = count[,keep]
pd_sub = pd_sub %>% select(cell_id, celltype) %>% left_join(gRNA_assign, by = "cell_id") %>% as.data.frame()
rownames(pd_sub) = as.vector(pd_sub$cell_id)
colnames(pd_sub) = c("cell_id", "cell_group", "embryo")

cds = new_cell_data_set(count_sub,
                        cell_metadata = pd_sub)

ccs = new_cell_count_set(cds, 
                         sample_group = "embryo", 
                         cell_group = "cell_group")

source("/help_code/power_analysis_hooke.R")
library(tibble)

effect_size_list = c(-0.75,-0.5,-0.25,-0.1,0.1,0.25,0.5,0.75)
num_embryos_list = c(10, 20, 30, 40, 50, 100, 200)
sample_size_list = c(50, 100, 200, 500, 1000)
success_ratio_list = c(0.5, 0.75, 1)

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1]) ### 1:10

for(sample_size in sample_size_list){
    result = NULL
    for(celltype_i in types){
        for(effect_i in effect_size_list){
            for(num_i in num_embryos_list){
                for(success_ratio_i in success_ratio_list){
                    print(paste0(sample_size, "/", celltype_i, "/", effect_i, "/", num_i))
                    res = run_comparison(ccs, 
                                         embryo_size = sample_size, 
                                         effect_size = (-1) * effect_i,
                                         num_embryos = num_i,
                                         cell_types = c(celltype_i),
                                         success_ratio = success_ratio_i,
                                         random.seed = kk, 
                                         num_threads = 4) 
                    result = rbind(result, data.frame(celltype = celltype_i,
                                                      effect_size = effect_i,
                                                      num_embryos = num_i,
                                                      success_ratio = success_ratio_i,
                                                      delta_log_abund = res$delta_log_abund[res$cell_group == celltype_i],
                                                      delta_p_value = res$delta_p_value[res$cell_group == celltype_i],
                                                      delta_q_value = res$delta_q_value[res$cell_group == celltype_i]))
                }
            }
        }
    }
    saveRDS(result, paste0(work_path, "/analysis/revision/power_analysis/iter_", kk, "/result_hooke_simulation_", sample_size, ".rds"))
}


###################################
### evaluating how variance of cell type compositions across samples will affect power

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1]) ### 1:10

rows = 22
mean_list = c(1000, 300, 100, 50)
sd_list = c(0.05, 0.25, 0.5)
sd_list_2 = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.25, 0.5, 1, 1.5, 2)

num_i = 30
effect_i = -0.5
success_ratio_i = 0.6
sample_size = 500

mat_all = NULL
for(mean_i in mean_list){
    for(sd_i in sd_list){
        mat_i = rnorm(rows, mean = mean_i, sd = mean_i * sd_i)
        mat_all = cbind(mat_all, mat_i)
    }
}
mat_orig = round(mat_all)

result = NULL
for(cnt in 1:length(sd_list_2)){
    mat_all = mat_orig
    mat_all[,5] = rnorm(rows, mean = mean_list[2], sd = mean_list[2] * sd_list_2[cnt])
    mat_all[mat_all < 0] = 1
    mat_all = t(mat_all)
    rownames(mat_all) = paste0("celltype_", 1:nrow(mat_all))
    colnames(mat_all) = paste0("sample_", 1:ncol(mat_all))
    ccs_new = monocle3::new_cell_data_set(mat_all)

    for(celltype_i in paste0("celltype_", 1:nrow(mat_all))){
        res = run_comparison(ccs_new, 
                             embryo_size = sample_size, 
                             effect_size = (-1) * effect_i,
                             num_embryos = num_i,
                             cell_types = c(celltype_i),
                             success_ratio = success_ratio_i,
                             random.seed = kk, 
                             num_threads = 4) 
        result = rbind(result, data.frame(sim = cnt,
                                          celltype = celltype_i,
                                          effect_size = effect_i,
                                          num_embryos = num_i,
                                          success_ratio = success_ratio_i,
                                          delta_log_abund = res$delta_log_abund[res$cell_group == celltype_i],
                                          delta_p_value = res$delta_p_value[res$cell_group == celltype_i],
                                          delta_q_value = res$delta_q_value[res$cell_group == celltype_i]))
    }
}

saveRDS(result, paste0(work_path, "/analysis/revision/power_analysis/iter_", kk, "/result_evaluate_variance.rds"))




#########################################
### Plotting the result by heatmap

library("gplots")
library(RColorBrewer)
library(viridis)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

effect_size_list = c(0.1,0.25,0.5,0.75)
num_embryos_list = c(10, 20, 30, 40, 50, 100, 200)
sample_size_list = c(50, 100, 200, 500, 1000)
success_ratio_list = c(0.5, 0.75, 1)

prob = readRDS(paste0(work_path, "/analysis/revision/power_analysis/celltype_prob.rds"))
colnames(prob) = c("celltype", "celltype_num", "celltype_prop")

dat = NULL
for(iter in 1:10){
    for(sample_size in sample_size_list){
        result = readRDS(paste0(work_path, "/analysis/revision/power_analysis/iter_", iter, "/result_hooke_simulation_", sample_size, ".rds"))
        result = result[result$effect_size < 0 & result$delta_log_abund < 0,]
        result$sample_size = sample_size
        dat = rbind(dat, result)
    }
}

dat = dat %>% left_join(prob[,c(1,3)], by = "celltype") %>% as.data.frame()
dat$log10_pval = -log10(dat$delta_p_value)
dat$effect_size = (-1)*dat$effect_size
dat$celltype = factor(dat$celltype, levels = as.vector(prob$celltype))

affect_factor_list = c("effect_size", "num_embryos", "success_ratio", "sample_size")

for (affect_factor in affect_factor_list) {
    message("Processing: ", affect_factor)
    df = dat %>% 
        group_by(celltype, !!sym(affect_factor)) %>% 
        summarize(mean_log10_pval = mean(log10_pval, na.rm = TRUE), .groups = "drop") %>%
        arrange(!!sym(affect_factor)) %>%
        dcast(celltype ~ get(affect_factor), value.var = "mean_log10_pval")
    rownames(df) = df$celltype
    mat = as.matrix(df[, -1, drop = FALSE])
    pdf(paste0(work_path, "/hooke_", affect_factor, ".pdf"), width = 8, height = 5)
    heatmap.2(
        mat,
        col = Colors,
        scale = "none",
        Rowv = FALSE,
        Colv = FALSE,
        key = TRUE,
        density.info = "none",
        trace = "none",
        cexRow = 0.5,
        cexCol = 0.5,
        margins = c(10, 5)
    )
    dev.off()
}


###############################################################
### Plotting the result by heatmap (only plotting the baseline)

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf"

library("gplots")
library(RColorBrewer)
library(viridis)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

effect_size_list = c(0.1,0.25,0.5,0.75)
num_embryos_list = c(10, 20, 30, 40, 50, 100, 200)
sample_size_list = c(50, 100, 200, 500, 1000)
success_ratio_list = c(0.5, 0.75, 1)

prob = readRDS(paste0(work_path, "/analysis/revision/power_analysis/celltype_prob.rds"))
colnames(prob) = c("celltype", "celltype_num", "celltype_prop")

dat = NULL
for(iter in 1:10){
    for(sample_size in sample_size_list){
        result = readRDS(paste0(work_path, "/analysis/revision/power_analysis/iter_", iter, "/result_hooke_simulation_", sample_size, ".rds"))
        result$delta_p_value[result$effect_size * result$delta_log_abund < 0] = 1
        result$sample_size = sample_size
        dat = rbind(dat, result)
    }
}

dat = dat %>% left_join(prob[,c(1,3)], by = "celltype") %>% as.data.frame()
dat$log10_pval = -log10(dat$delta_p_value)
dat$effect_size = abs(dat$effect_size)
dat$celltype = factor(dat$celltype, levels = as.vector(prob$celltype))
dat = dat[,c("celltype", "effect_size", "num_embryos", "success_ratio", "sample_size", "log10_pval")]

affect_factor_list = c("effect_size", "num_embryos", "success_ratio", "sample_size")

for (affect_factor in affect_factor_list) {
    message("Processing: ", affect_factor)
    if(affect_factor == "effect_size") {
        df = dat %>% filter(num_embryos == 30, success_ratio == 0.75, sample_size == 500)
    } else if (affect_factor == "num_embryos") {
        df = dat %>% filter(effect_size == 0.5, success_ratio == 0.75, sample_size == 500)
    } else if (affect_factor == "success_ratio") {
        df = dat %>% filter(effect_size == 0.5, num_embryos == 30, sample_size == 500)
    } else {
        df = dat %>% filter(effect_size == 0.5, num_embryos == 30, success_ratio == 0.75)
    }
    df = df %>% 
        group_by(celltype, !!sym(affect_factor)) %>% 
        summarize(mean_log10_pval = mean(log10_pval, na.rm = TRUE), .groups = "drop") %>%
        arrange(!!sym(affect_factor)) %>%
        dcast(celltype ~ get(affect_factor), value.var = "mean_log10_pval")
    rownames(df) = df$celltype
    mat = as.matrix(df[, -1, drop = FALSE])
    pdf(paste0(work_path, "/hooke_", affect_factor, ".pdf"), width = 8, height = 5)
    heatmap.2(
        mat,
        col = Colors,
        scale = "none",
        Rowv = FALSE,
        Colv = FALSE,
        key = TRUE,
        density.info = "none",
        trace = "none",
        cexRow = 0.5,
        cexCol = 0.5,
        margins = c(10, 5)
    )
    dev.off()
}




#########################################
### Plotting the result of evaluating variance of cell type compositions across samples


library("gplots")
library(RColorBrewer)
library(viridis)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

mean_list = c(1000, 300, 100, 50)
sd_list = c(0.05, 0.25, 0.5)
sd_list_2 = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.25, 0.5, 1, 1.5, 2)

num_i = 30
effect_i = -0.5
success_ratio_i = 0.6
sample_size = 500

dat = NULL
for(iter in 1:10){
    dat_i = readRDS(paste0(work_path, "/analysis/revision/power_analysis/iter_", iter, "/result_evaluate_variance.rds"))
    dat = rbind(dat, dat_i)
}
dat$log10_pval = -log10(dat$delta_p_value)

mean_list = c(1000, 300, 100, 50)
sd_list = c(0.05, 0.25, 0.5)
celltype_list = data.frame(celltype = paste0("celltype_", 1:12),
    mean_num = rep(mean_list, each = 3),
    sd_by_mean = rep(sd_list, 4))

df = dat[dat$celltype == "celltype_5" & dat$delta_log_abund < 0,] %>%
     group_by(sim) %>% summarize(mean_log10_pval = mean(log10_pval))
df$sd_by_mean = factor(sd_list_2, levels = sd_list_2)

p = ggplot() +
    geom_point(data = df, aes(x = sd_by_mean, y = mean_log10_pval)) +
    geom_line(data = df, aes(x = sd_by_mean, y = mean_log10_pval, group = 1)) +
    labs(x="SD_by_Mean", y="Mean_log10_Pval", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0(work_path, "/result_evaluate_variance.pdf"),p,width = 5, height = 3)








