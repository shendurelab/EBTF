
##################################
### Investigating the mosaic issue
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tf
### validated_125TF.cell_gene.rds

#################################
### Step-1: cell number comparing

pd = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned_pd.rds"))
pd$cell = gsub("EB21_CRISPRcut_125TF_GEx_", "", as.vector(pd$cell_id))
cell_gene = readRDS("validated_125TF.cell_gene.rds")

pd$replicate = rep("rep1", nrow(pd))
pd$replicate[pd$experiment_id %in% paste0("EB21_CRISPRcut_125TF_GEx_EB_", c(17:32))] = "rep2"
cell_gene = cell_gene %>% left_join(pd[,c("cell","celltype","replicate")], by = "cell")

df = cell_gene %>% group_by(barcode, replicate) %>% tally() %>%
    dcast(barcode~replicate, fill = 0) %>% mutate(log2_rep1 = log2(rep1+1), log2_rep2 = log2(rep2+1))

cor.test(df$log2_rep1, df$log2_rep2)

p = ggplot() + 
    geom_point(data = df, aes(x = log2_rep1, y = log2_rep2)) + 
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/mosaic_plot/gRNA_number_two_rep2.pdf"))
print(p)
dev.off()


##########
### Step-2 

dat = read.table("/net/shendure/vol10/projects/silvia/EBs/scaled_125TFs/nobackup/clonal/gRNA_counts.txt", header=T, as.is=T)

df = dat %>% select(STD_10k_1, STD_10k_2) %>% filter(STD_10k_1 > 1000, STD_10k_2 > 1000) %>%
    mutate(log2_rep1 = log2(STD_10k_1), log2_rep2 = log2(STD_10k_2))

cor.test(df$log2_rep1, df$log2_rep2)

p = ggplot() + 
    geom_point(data = df, aes(x = log2_rep1, y = log2_rep2)) + 
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/mosaic_plot/gRNA_number_two_ESs.pdf"))
print(p)
dev.off()

res = NULL
for(i in c("100", "1k", "5k", "10k", "50k", "100k")){
    dat$rep1 = dat[[paste0("STD_", i, "_1")]]
    dat$rep2 = dat[[paste0("STD_", i, "_2")]]
    
    df = dat %>% select(rep1, rep2) %>% filter(rep1 > 1000, rep2 > 1000) %>%
        mutate(log2_rep1 = log2(rep1), log2_rep2 = log2(rep2))
    
    x = cor.test(df$log2_rep1, df$log2_rep2, method = "spearman")
    res = rbind(res, data.frame(cell_num = i,
                                estimate = x$estimate,
                                pval = x$p.value))
    
    p = ggplot() + 
        geom_point(data = df, aes(x = log2_rep1, y = log2_rep2)) + 
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
    pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/mosaic_plot/gRNA_number_two_ESs_", i,".pdf"))
    print(p)
    dev.off()
}


df = dat %>% select(EB_1, EB_2) %>% filter(EB_1 > 1000, EB_2 > 1000) %>%
    mutate(log2_EB1 = log2(EB_1), log2_EB2 = log2(EB_2))

cor.test(df$log2_EB1, df$log2_EB2)

p = ggplot() + 
    geom_point(data = df, aes(x = log2_EB1, y = log2_EB2)) + 
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/mosaic_plot/gRNA_number_EB1_EB2.pdf"))
print(p)
dev.off()




# Stacked + percent
df = data.frame(num = c(as.matrix(dat[,c(paste0("EB_", 1:12), "STD_10k_1", "STD_10k_2")])),
                group = rep(c(paste0("EB_", 1:12), "STD_10k_1", "STD_10k_2"), each = nrow(dat)),
                guide = rep(dat$gRNA, times = 14))

df$group = factor(df$group, levels = c(paste0("EB_", 1:12), "STD_10k_1", "STD_10k_2"))

p = df %>% 
    #filter(guide %in% unique(dat$gRNA)[1:10]) %>%
    group_by(group, guide) %>%
    summarize(n = sum(num)) %>%
    ggplot(aes(fill=guide, y=n, x=group)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    labs(x="", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/mosaic_plot/stacked.pdf"),
           dpi = 300,
           height  = 5, 
           width = 12)


##########
### Step-3

# for each sample, the distribution of read numbers from the top1 guide, top2 guide, ...

sample_list = c(paste0("EB_", 1:12), "STD_10k_1", "STD_10k_2")

df = NULL

for(sample_i in sample_list){
    x = dat[,sample_i]
    x = sort(x, decreasing=T)
    x = data.frame(num = x,
                   cum_num = cumsum(x),
                   guide = c(1:456),
                   group = sample_i)
    df = rbind(df, x)
}

df_sum = df %>% group_by(group) %>% summarize(num_sum = sum(num))
df = df %>% 
    left_join(df_sum, by = "group") %>%
    mutate(frac = 100*cum_num/num_sum) 

col_value = c(rep("grey50", 12), rep("red", 2))

ggplot(df, aes(x = guide, y = frac, color = group)) + 
    #geom_point() + 
    geom_line(linewidth = 1) +
    labs(title = "", x = "guides", y = "cumulative % of reads") +
    theme_classic(base_size = 15) +
    #scale_x_continuous(limits = c(1, 456), breaks = c(1:456)) +
    theme(legend.position="none") +
    scale_color_manual(values=col_value) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/mosaic_plot/line_plot.pdf"),
           dpi = 300,
           height  = 5, 
           width = 8)









