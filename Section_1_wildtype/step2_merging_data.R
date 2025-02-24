
#####################################################
### Merging three timepoints and performing embedding
### Chengxiang Qiu
### Feb-20, 2025

###############################################################
### Step-1: merging data from 3 timepoints (day7, day14, day21)

source("help_script.R")
mouse_gene <- read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

run_id = "day7"
obj_1 = readRDS(paste0(paste0(work_path, "/processing/mEB_time_course_GEx/", run_id, "/obj_", "mEB_time_course_GEx_", run_id, ".rds")))
run_id = "day14"
obj_2 = readRDS(paste0(paste0(work_path, "/processing/mEB_time_course_GEx/", run_id, "/obj_", "mEB_time_course_GEx_", run_id, ".rds")))
run_id = "day21"
obj_3 = readRDS(paste0(paste0(work_path, "/processing/mEB_time_course_GEx/", run_id, "/obj_", "mEB_time_course_GEx_", run_id, ".rds")))

obj = merge(obj_1, c(obj_2, obj_3))

count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]
obj = CreateSeuratObject(count, meta.data = pd)

saveRDS(count, paste0(work_path, "/analysis/mEB_time_course_GEx/count.rds"))
saveRDS(pd, paste0(work_path, "/analysis/mEB_time_course_GEx/pd.rds"))

##################################
### Step-2: performing sctransform

count = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/pd.rds"))

obj = CreateSeuratObject(count, meta.data = pd)
obj = SCTransform(object = obj, verbose = FALSE)
obj = RunPCA(object = obj, npcs = 30, verbose = FALSE)
obj = RunUMAP(object = obj, dims = 1:30, min.dist = 0.3, verbose = FALSE)
obj = FindNeighbors(object = obj, dims = 1:30, verbose = FALSE, reduction = "pca")
obj = FindClusters(object = obj, resolution = 1, verbose = FALSE)
obj = FindClusters(object = obj, resolution = 2, verbose = FALSE)
saveRDS(obj, paste0(work_path, "/analysis/mEB_time_course_GEx/obj_sctransform.rds"))

pd = data.frame(obj[[]])
pd$UMAP_1 = as.vector(Embeddings(obj, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/analysis/mEB_time_course_GEx/obj_sctransform_pd.rds"))

for(i in unique(pd$experiment_id)){
    print(i)
    try(ggplot(pd) +
            geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.1, color = "grey80") +
            geom_point(data = subset(pd, experiment_id == i),
                       aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.15) +
            geom_point(data = subset(pd, experiment_id == i),
                       aes(x = UMAP_1, y = UMAP_2, color = experiment_id), size=0.1) +
            theme_void() +
            theme(legend.position="none") +
            scale_color_manual(values=mEB_time_course_color_code) +
            ggsave(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/obj_sctransform_", i, ".png"),
                   dpi = 300,
                   height  = 5, 
                   width = 5), silent = TRUE)
}

anno = rep(NA, nrow(pd))
anno[pd$SCT_snn_res.1 %in% c(4)] = "Primordial germ cells"
anno[pd$SCT_snn_res.1 %in% c(0,6,10,13)] = "Epiblast and primitive streak"
anno[pd$SCT_snn_res.1 %in% c(1)] = "Neuroectoderm"
anno[pd$SCT_snn_res.1 %in% c(7)] = "Early neurons"
anno[pd$SCT_snn_res.1 %in% c(9)] = "Surface ectoderm"
anno[pd$SCT_snn_res.1 %in% c(12,17)] = "Definitive endoderm"
anno[pd$SCT_snn_res.1 %in% c(3,5,15,16,19)] = "Extraembryonic visceral endoderm"
anno[pd$SCT_snn_res.1 %in% c(2)] = "Parietal endoderm"
anno[pd$SCT_snn_res.1 %in% c(11)] = "Nascent mesoderm"
anno[pd$SCT_snn_res.1 %in% c(8,14)] = "Lateral plate mesoderm"
anno[pd$SCT_snn_res.1 %in% c(21)] = "Cardiomyocytes"
anno[pd$SCT_snn_res.1 %in% c(20)] = "Endothelial cells"
anno[pd$SCT_snn_res.1 %in% c(18)] = "Blood progenitors"

anno[pd$RNA_snn_res.1 %in% c(17)] = "Notochord"
anno[pd$RNA_snn_res.1 %in% c(20)] = "Erythroid cells"

pd$celltype = as.vector(anno)
saveRDS(pd, paste0(work_path, "/analysis/mEB_time_course_GEx/obj_sctransform_pd.rds"))


try(ggplot(pd) +
        geom_point(aes(x = UMAP_1, y = UMAP_2), size=0.3, color = "black") +
        geom_point(aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.2) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/mEB_time_course_GEx_sctransform_celltype.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

######################################################
### Step-3: cell compositions change across timepoints

EB_celltype_color_code = EB_celltype_color_code[names(EB_celltype_color_code) %in% pd$celltype]
pd$day = gsub("mEB_time_course_GEx_", "", as.vector(pd$experiment_id))
pd$day = factor(pd$day, levels = paste0("day", c(7, 14, 21)))
pd$celltype = factor(pd$celltype, levels = names(EB_celltype_color_code))

# Stacked + percent
p = pd %>% 
    group_by(day, celltype) %>%
    tally() %>%
    ggplot(aes(fill=celltype, y=n, x=day)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=EB_celltype_color_code) +
    labs(x="", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ggsave(paste0(work_path, "/analysis/mEB_time_course_GEx/plot/cell_composition.pdf"),
           dpi = 300,
           height  = 10, 
           width = 5)


celltype_regroup = rep(NA, nrow(pd))
celltype_regroup[pd$celltype %in% c("Primordial germ cells",
                                    "Epiblast and primitive streak")] = "Epiblast"
celltype_regroup[pd$celltype %in% c("Neuroectoderm",
                                    "Early neurons")] = "Neuroectoderm"
celltype_regroup[pd$celltype %in% c("Notochord")] = "Notochord"
celltype_regroup[pd$celltype %in% c("Surface ectoderm")] = "Surface ectoderm"
celltype_regroup[pd$celltype %in% c("Definitive endoderm")] = "Def. endoderm"
celltype_regroup[pd$celltype %in% c("Extraembryonic visceral endoderm",
                                    "Parietal endoderm")] = "ExE endoderm"
celltype_regroup[pd$celltype %in% c("Nascent mesoderm",
                                    "Lateral plate mesoderm")] = "Mesoderm"
celltype_regroup[pd$celltype %in% c("Nascent mesoderm",
                                    "Lateral plate mesoderm",
                                    "Cardiomyocytes")] = "Mesoderm"
celltype_regroup[pd$celltype %in% c("Endothelial cells",
                                    "Blood progenitors",
                                    "Erythroid cells")] = "Blood"
pd$celltype_regroup = as.vector(celltype_regroup)
saveRDS(pd, paste0(work_path, "/analysis/mEB_time_course_GEx/obj_sctransform_pd.rds"))

pd_x = pd %>% group_by(celltype_regroup, experiment_id) %>% tally() %>%
    left_join(pd %>% group_by(experiment_id) %>% tally() %>% rename(total_n = n), by = "experiment_id") %>%
    mutate(pct = 100*n/total_n) %>% as.data.frame()


