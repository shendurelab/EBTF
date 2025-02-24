
#####################################################################
### Merging data from pilot experiment to perform dimension reduction
### Chengxiang Qiu
### Feb-20, 2025

##############################################################################
### Step-1: merging data from wildtype-day21, six batches of pilot experiments

### mEB_time_course_GEx_day21
### rep1_EB21_CRISPRcut_sorted *** excluded due to low number of cells (n = 2,699) and >60% of sgRNAs were assigned with a specific NTC barcode.
### rep1_EB21_CRISPRcut_unsorted
### rep1_EB21_CRISPRi_sorted
### rep1_EB21_CRISPRi_unsorted
### rep2_EB21_CRISPRcut_sorted
### rep2_EB21_CRISPRi_sorted

source("help_script.R")
mouse_gene <- read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

run_id = "day21"
obj_1 = readRDS(paste0(paste0(work_path, "/processing/mEB_time_course_GEx/", run_id, "/obj_", "mEB_time_course_GEx_", run_id, ".rds")))
run_id = "rep1_EB21_CRISPRcut_unsorted"
obj_2 = readRDS(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/obj_", "pilot_screen_GEx_", run_id, ".rds"))
run_id = "rep1_EB21_CRISPRi_unsorted"
obj_3 = readRDS(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/obj_", "pilot_screen_GEx_", run_id, ".rds"))
run_id = "rep2_EB21_CRISPRcut_sorted"
obj_4 = readRDS(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/obj_", "pilot_screen_GEx_", run_id, ".rds"))
run_id = "rep2_EB21_CRISPRi_sorted"
obj_5 = readRDS(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/obj_", "pilot_screen_GEx_", run_id, ".rds"))
run_id = "rep1_EB21_CRISPRi_sorted"
obj_6 = readRDS(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/obj_", "pilot_screen_GEx_", run_id, ".rds"))

obj = merge(obj_1, c(obj_2, obj_3, obj_4, obj_5, obj_6))

count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])

mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "M")) &
                                mouse_gene$gene_type %in% c("protein_coding", "lincRNA"),]
count = count[rownames(count) %in% as.vector(mouse_gene_sub$gene_ID),]
obj = CreateSeuratObject(count, meta.data = pd)

saveRDS(count, paste0(work_path, "/analysis/pilot_screen_GEx/count.rds"))
saveRDS(pd, paste0(work_path, "/analysis/pilot_screen_GEx/pd.rds"))


################################################
### Step-2: align by correcting the batch effect


count = readRDS(paste0(work_path, "/analysis/pilot_screen_GEx/count.rds"))
pd = readRDS(paste0(work_path, "/analysis/pilot_screen_GEx/pd.rds"))
obj = CreateSeuratObject(count, meta.data = pd)

group = rep("mEB_time_course_GEx_day21", nrow(pd))
group[pd$experiment_id %in% c("pilot_screen_GEx_rep1_EB21_CRISPRcut_unsorted", "pilot_screen_GEx_rep1_EB21_CRISPRi_unsorted", "pilot_screen_GEx_rep1_EB21_CRISPRi_sorted")] = "pilot_screen_GEx_rep1"
group[pd$experiment_id %in% c("pilot_screen_GEx_rep2_EB21_CRISPRcut_sorted", "pilot_screen_GEx_rep2_EB21_CRISPRi_sorted")] = "pilot_screen_GEx_rep2"
print(table(group))

obj$group = as.vector(group)
obj_processed = doClusterSeurat(obj)
obj_processed = FindClusters(object = obj_processed, resolution = 2, verbose = FALSE)
saveRDS(obj_processed, paste0(work_path, "/analysis/pilot_screen_GEx/obj_aligned.rds"))

pd = data.frame(obj_processed[[]])
pd$UMAP_1 = as.vector(Embeddings(obj_processed, reduction = "umap")[,1])
pd$UMAP_2 = as.vector(Embeddings(obj_processed, reduction = "umap")[,2])
pd$group = NULL
saveRDS(pd, paste0(work_path, "/analysis/pilot_screen_GEx/obj_aligned_pd.rds"))


###################################################
### Step-3: making plot based on integration result

pd = readRDS(paste0(work_path, "/analysis/pilot_screen_GEx/obj_aligned_pd.rds"))

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
            scale_color_manual(values=pilot_screen_GEx_color_code) +
            ggsave(paste0(work_path, "/analysis/pilot_screen_GEx/plot/obj_aligned_", i, ".png"),
                   dpi = 300,
                   height  = 5, 
                   width = 5), silent = TRUE)
}

obj_processed = readRDS(paste0(work_path, "/analysis/pilot_screen_GEx/obj_aligned.rds"))
pca_coor = Embeddings(obj_processed, reduction = "pca")

pd_1 = pd[pd$experiment_id != "mEB_time_course_GEx_day21",]
pd_2 = pd[pd$experiment_id == "mEB_time_course_GEx_day21",]
pca_coor_1 = pca_coor[pd$experiment_id != "mEB_time_course_GEx_day21",]
pca_coor_2 = pca_coor[pd$experiment_id == "mEB_time_course_GEx_day21",]

k.param = 10; nn.method = "rann"; nn.eps = 0; annoy.metric = "euclidean"
nn.ranked = Seurat:::NNHelper(
    data = pca_coor_2,
    query = pca_coor_1,
    k = k.param,
    method = nn.method,
    searchtype = "standard",
    eps = nn.eps,
    metric = annoy.metric)
nn.ranked = Indices(object = nn.ranked)
nn_matrix = nn.ranked

pd_x = readRDS(paste0(work_path, "/analysis/mEB_time_course_GEx/obj_sctransform_pd.rds"))
pd_x = pd_x[row.names(pd_2),]
pd_2$celltype = as.vector(pd_x$celltype)

resultA = NULL
for(i in 1:k.param){
    print(i)
    resultA = cbind(resultA, as.vector(pd_2$celltype)[as.vector(nn_matrix[,i])])
}
resultB = lapply(1:nrow(resultA), function(x){
    return(names(sort(table(resultA[x,]), decreasing = T))[1])
})
pd_1$celltype = as.vector(unlist(resultB))

pd = rbind(pd_1, pd_2)
pd = pd[colnames(obj_processed),]
saveRDS(pd, paste0(work_path, "/analysis/pilot_screen_GEx/obj_aligned_pd.rds"))

try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), size=0.2, color = "grey90") +
        geom_point(data = subset(pd, experiment_id == "mEB_time_course_GEx_day21"),
                   aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.6) +
        geom_point(data = subset(pd, experiment_id == "mEB_time_course_GEx_day21"),
                   aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.4) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/pilot_screen_GEx/plot/pilot_screen_GEx_EB21_wt_celltype.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)

try(ggplot() +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2), color = "black", size=0.2) +
        geom_point(data = pd, aes(x = UMAP_1, y = UMAP_2, color = celltype), size=0.15) +
        theme_void() +
        theme(legend.position="none") +
        scale_color_manual(values=EB_celltype_color_code) +
        ggsave(paste0(work_path, "/analysis/pilot_screen_GEx/plot/pilot_screen_celltype.png"),
               dpi = 300,
               height  = 5, 
               width = 5), silent = TRUE)


####################################################
### Step-4: cell compositions change across datasets

EB_celltype_color_code = EB_celltype_color_code[names(EB_celltype_color_code) %in% pd$celltype]
pd$celltype = factor(pd$celltype, levels = names(EB_celltype_color_code))

pd$experiment_id = factor(pd$experiment_id, levels = c("mEB_time_course_GEx_day21",
                                                       "pilot_screen_GEx_rep1_EB21_CRISPRcut_unsorted",
                                                       "pilot_screen_GEx_rep1_EB21_CRISPRi_unsorted",
                                                       "pilot_screen_GEx_rep2_EB21_CRISPRcut_sorted",
                                                       "pilot_screen_GEx_rep1_EB21_CRISPRi_sorted",
                                                       "pilot_screen_GEx_rep2_EB21_CRISPRi_sorted"))

# Stacked + percent
p = pd %>% 
    group_by(experiment_id, celltype) %>%
    tally() %>%
    ggplot(aes(fill=celltype, y=n, x=experiment_id)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=EB_celltype_color_code) +
    labs(x="", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggsave(paste0(work_path, "/analysis/pilot_screen_GEx/plot/cell_composition.pdf"),
           dpi = 300,
           height  = 12, 
           width = 7)

### cell compositions change across datasets (only focusing on NTCs)
### looking at if the high % of Epiblast in the perturbation dataset is due to sorting bias

pd = readRDS(paste0(work_path, "/analysis/pilot_screen_GEx/obj_aligned_pd.rds"))

run_list = c("rep2_EB21_CRISPRcut_sorted_gRNA", 
             "rep1_EB21_CRISPRi_sorted_gRNA", 
             "rep2_EB21_CRISPRi_sorted_gRNA")

NTC_list = NULL
for(run_id in run_list){
    dat = read.table(paste0(work_path, "/processing/pilot_screen_sgRNA/", run_id, "/cell_gene.txt"), header=T, as.is=T, sep="\t")
    dat$cell = paste0("pilot_screen_GEx_", as.vector(dat$cell))
    print(nrow(dat))
    print(sum(dat$cell %in% rownames(pd)))
    
    x = as.vector(dat$gene)
    x[dat$gene %in% c("NTC.control-random-chr3-116819213-116820756-1",
                      "NTC.control-random-chr3-116819213-116820756-2",
                      "NTC.control-random-chr5-31878254-31878457-3",
                      "NTC.control-random-chr5-31878254-31878457-4",
                      "NTC.control-random-chr7-114416771-114416984-5",
                      "NTC.control-random-chr7-114416771-114416984-6",
                      "NTC.CRISPRCUT-NONTARGETING",
                      "NTC.CRISPRI-NONTARGETING")] = "NTC"
    dat$gene = as.vector(x)    
    
    NTC_list = c(NTC_list, as.vector(dat$cell[dat$gene == "NTC"]))
}

pd_sub = pd[pd$cell_id %in% NTC_list,]
pd = pd_sub

EB_celltype_color_code = EB_celltype_color_code[names(EB_celltype_color_code) %in% pd$celltype]
pd$celltype = factor(pd$celltype, levels = names(EB_celltype_color_code))

pd$experiment_id = factor(pd$experiment_id, levels = c("pilot_screen_GEx_rep2_EB21_CRISPRcut_sorted",
                                                       "pilot_screen_GEx_rep1_EB21_CRISPRi_sorted",
                                                       "pilot_screen_GEx_rep2_EB21_CRISPRi_sorted"))

# Stacked + percent
p = pd %>% 
    group_by(experiment_id, celltype) %>%
    tally() %>%
    ggplot(aes(fill=celltype, y=n, x=experiment_id)) + 
    geom_bar(position="fill", stat="identity", width = 0.8) +
    scale_fill_manual(values=EB_celltype_color_code) +
    labs(x="", y="% of cells", title="") +
    theme_classic(base_size = 15) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggsave(paste0(work_path, "/analysis/pilot_screen_GEx/plot/cell_composition_NTC.pdf"),
           dpi = 300,
           height  = 12, 
           width = 7)







