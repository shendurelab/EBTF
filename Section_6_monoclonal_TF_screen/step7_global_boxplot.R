

print("Loading packages for regular data analysis, e.g. dplyr")
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))

print("Loading packages for plotting, e.g. ggplot2")
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(gridExtra))
suppressMessages(library(viridis))
suppressMessages(library(gplots))
suppressMessages(library(patchwork))

work_path = ""

EB_celltype_color_code = c("Primordial germ cells" = "#9c476b",
                           "ESCs (2-cell state)" = "#379884",
                           "Epiblast" = "#5dbb4d",
                           "Epiblast and primitive streak" = "#5dbb4d",
                           "Primitive streak" = "#d3756f",
                           "Neuroectoderm" = "#7459c6",
                           "Surface ectoderm" = "#c845a1",
                           "Notochord" = "#46873f",
                           "Endoderm" = "#bc76d5",
                           "Definitive endoderm" = "#bc76d5",
                           "Extraembryonic visceral endoderm" = "#d59a35",
                           "Parietal endoderm" = "#6c80c9",
                           "Paraxial mesoderm" = "#dc87ba",
                           "Nascent mesoderm" = "#d05130",
                           "Lateral plate mesoderm" = "#4eacd7",
                           "Cardiomyocytes" = "#73732a",
                           "Endothelial cells" = "#5dc598",
                           "Blood progenitors" = "#d681b2",
                           "Erythroid cells" = "#9cb26a",
                           "Early neurons" = "#a8b439",
                           "Floor plate" = "#9f5f2b",
                           "Eye field" = "#ac4c55")

############
### Arraryed

dat = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep1/hooke_ccs_pd.rds"))
pd = dat[['pd']]
cell_freq = dat[['cell_freq']]

cell_freq = 100*t(t(cell_freq)/colSums(cell_freq))

df = melt(as.matrix(cell_freq))
colnames(df) = c("celltype", "sample", "freq")
df = df %>% left_join(pd[,c('sample','target')], by = 'sample')

gene_list = c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1","NONTARGETING")
celltype_list = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% df$celltype]

df$target = factor(df$target, levels = gene_list)
df$celltype = factor(df$celltype, levels = celltype_list)

p = ggplot(df, aes(target, freq, fill = celltype)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5) + 
    labs(x="", y="% of cells") +
    facet_wrap(~celltype, ncol = 1) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_fill_manual(values=EB_celltype_color_code) +
    facet_wrap(~celltype, ncol = 1, scales = "free_y") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 45, hjust = 1), axis.text.y = element_text(color="black"))
ggsave("~/share/big_boxplot_rep1.pdf", p, height = 12, width = 5)




##########
### Pooled

dat = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/celltype_changes/rep2/hooke_ccs_pd.rds"))
pd = dat[['pd']]
cell_freq = dat[['cell_freq']]

cell_freq = 100*t(t(cell_freq)/colSums(cell_freq))

df = melt(as.matrix(cell_freq))
colnames(df) = c("celltype", "sample", "freq")
df = df %>% left_join(pd[,c('sample','target')], by = 'sample')

gene_list = c("T","Lhx1","Hand1","Gata6","Hes5","Six3","Lmo2","Hhex","Carm1","NONTARGETING")
celltype_list = names(EB_celltype_color_code)[names(EB_celltype_color_code) %in% df$celltype]

df$target = factor(df$target, levels = gene_list)
df$celltype = factor(df$celltype, levels = celltype_list)

p = ggplot(df, aes(target, freq, fill = celltype)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5) + 
    labs(x="", y="% of cells") +
    facet_wrap(~celltype, ncol = 1) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    scale_fill_manual(values=EB_celltype_color_code) +
    facet_wrap(~celltype, ncol = 1, scales = "free_y") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black", angle = 45, hjust = 1), axis.text.y = element_text(color="black"))
ggsave("~/share/big_boxplot_rep2.pdf", p, height = 12, width = 5)

