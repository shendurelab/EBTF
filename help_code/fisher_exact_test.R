
#####################################
### everytime, only do 1 TF (n = 820)

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])
target_i = target_cell_num$gene[kk]
print(target_i)

cell_gene_target = cell_gene %>% filter(gene == target_i) %>%
    select(cell, celltype) %>% unique() %>%
    group_by(celltype) %>% tally() %>% rename(target_n = n)

res = NULL
iter_time = 5000

for(cnt in c(1:iter_time)){
    
    print(cnt)
    
    cell_gene_NTC = cell_gene %>% filter(barcode %in% NTC_cell_num$barcode[sample(1:nrow(NTC_cell_num), 3)]) %>%
        select(cell, celltype) %>% unique() %>%
        group_by(celltype) %>% tally() %>% rename(NTC_n = n)
    
    compare_i = cell_gene_target %>% full_join(cell_gene_NTC, by = "celltype") %>%
        mutate_if(is.numeric,coalesce,0)
    
    for(j in 1:nrow(compare_i)){
        a = compare_i$target_n[j]
        b = compare_i$NTC_n[j]
        fit = fisher.test(matrix(c(a, b, sum(compare_i$target_n) - a, sum(compare_i$NTC_n) - b), 2, 2))
        res = rbind(res, data.frame(celltype = compare_i$celltype[j],
                                    odds_ratio = fit$estimate,
                                    p_val = fit$p.val))
    }
    
}

rownames(res) = NULL
saveRDS(res, paste0(work_path, "/analysis/all_TF_screen_GEx/fisher_test/", target_i, ".rds"))


#############################################################
### adjusting p-value across all the tests on different genes

source("~/work/scripts/utils.R")
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf"

res = NULL
file_list = dir(paste0(work_path, "/analysis/all_TF_screen_GEx/fisher_test/"))
for(cnt in 1:length(file_list)){
    print(paste0(cnt, "/", length(file_list)))
    res_i = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/fisher_test/", file_list[cnt]))
    res = c(res, as.vector(res_i$p_val))
}
res_adjust = p.adjust(res, method = "fdr")
saveRDS(res_adjust, paste0(work_path, "/analysis/all_TF_screen_GEx/fisher_test_fdr.rds"))

for(cnt in 1:length(file_list)){
    print(paste0(cnt, "/", length(file_list)))
    res_i = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/fisher_test/", file_list[cnt]))
    res_i$fdr = res_adjust[1:nrow(res_i)]
    res_adjust = res_adjust[(nrow(res_i) + 1):length(res_adjust)]
    saveRDS(res_i, paste0(work_path, "/analysis/all_TF_screen_GEx/fisher_test/", file_list[cnt]))
}



#####################################
### summarize the results

res_all = NULL

cnt = 1
for(target_i in as.vector(target_cell_num$gene)){
    print(paste0(cnt, "/", target_i))
    cnt = cnt + 1
    
    res = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/fisher_test/", target_i, ".rds"))
    
    res_i = res %>% group_by(celltype) %>% 
        summarize(median_or = median(odds_ratio), median_p = median(p_val), median_fdr = median(fdr)) %>% 
        arrange(desc(median_or)) %>% mutate(target = target_i)
    
    res_all = rbind(res_all, res_i)
    
}

res_all = res_all %>% arrange(median_fdr)
res_all$log2_median_or = log2(res_all$median_or)
saveRDS(res_all, paste0(work_path, "/analysis/all_TF_screen_GEx/fisher_test_result.rds"))

res_sig = res_all %>% filter(median_fdr < 0.01 & abs(log2_median_or) > 1)
saveRDS(res_sig, paste0(work_path, "/analysis/all_TF_screen_GEx/fisher_test_result_sig.rds"))

res_all$log2_median_or[is.infinite(res_all$log2_median_or) & res_all$log2_median_or > 0] = 8
res_all$log2_median_or[is.infinite(res_all$log2_median_or) & res_all$log2_median_or < 0] = -8
res_all$log10_median_fdr = -log10(res_all$median_fdr)

ggplot() +
    geom_point(data = res_all, aes(x = log2_median_or, y = log10_median_fdr), size=0.1, color = "grey80") +
    geom_point(data = subset(res_all, median_fdr < 0.01 & abs(log2_median_or) > 1), aes(x = log2_median_or, y = log10_median_fdr, color = celltype), size=0.5) +
    theme_classic(base_size = 10) +
    labs(x="log2 median odds ratio", y="-log10 median FDR") +
    theme(legend.position="none") +
    geom_vline(xintercept = c(-1,1), linetype="dotted") +
    geom_hline(yintercept = 2, linetype="dotted") +
    scale_color_manual(values=EB_celltype_color_code) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
ggsave(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/fisher_test_vacano_plot.png"),
       dpi = 300,               
       height  = 5, 
       width = 5)

p = ggplot() +
    geom_point(data = res_all, aes(x = log2_median_or, y = log10_median_fdr), size=0.1, color = "grey80") +
    geom_point(data = subset(res_all, median_fdr < 0.01 & abs(log2_median_or) > 1), aes(x = log2_median_or, y = log10_median_fdr, color = celltype), size=0.5) +
    theme_classic(base_size = 10) +
    labs(x="log2 median odds ratio", y="-log10 median FDR") +
    theme(legend.position="none") +
    geom_vline(xintercept = c(-1,1), linetype="dotted") +
    geom_hline(yintercept = 2, linetype="dotted") +
    scale_color_manual(values=EB_celltype_color_code) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))
pdf(paste0(work_path, "/analysis/all_TF_screen_GEx/plot/fisher_test_vacano_plot.pdf"),5,5)
print(p)
dev.off()


### select some TFs for UMAP visualization
res_sub = res_sig %>% filter(log2_median_or > 0) %>% group_by(celltype) %>%
    slice_min(order_by = median_fdr, n = 1) %>% as.data.frame()

