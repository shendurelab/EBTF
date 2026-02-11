

#########################
### Dataset-1: big screen

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

obj = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned.rds"))
pd = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned_pd.rds"))
pd$cell = gsub("all_TF_screen_GEx_", "", as.vector(pd$cell_id))

obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1000)
exp = GetAssayData(obj, slot = "counts")[VariableFeatures(obj),]
emb = as.matrix(pd[,c("UMAP_1","UMAP_2")])
pd = pd[,c("cell","celltype","experiment_id")]
fd = mouse_gene[rownames(exp),]

cds = new_cell_data_set(exp,
                        cell_metadata = pd,
                        gene_metadata = fd)
cds = preprocess_cds(cds)
reducedDims(cds)$UMAP = emb
cds = cluster_cells(cds)
cds = learn_graph(cds)

saveRDS(cds, paste0(work_path, "/analysis/revision/pseudotime/cds_big_screen.rds"))

### START ###
### using monocle3 to estimate pseudotime (locally)
source("~/work/scripts/utils.R")
cds = readRDS("~/work/sam_tf/revision/pseudotime/cds_big_screen.rds")
cds = order_cells(cds)
pd = data.frame(pData(cds))
pd$pseudotime = cds@principal_graph_aux[["UMAP"]]$pseudotime
saveRDS(pd, "~/work/sam_tf/revision/pseudotime/cds_big_screen_pd.rds")

p = plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
ggsave("~/work/sam_tf/revision/pseudotime/cds_big_screen_pseudotime.png", dpi=300, width=7, height=5)
### END ###

pd = readRDS(paste0(work_path, "/analysis/revision/pseudotime/cds_big_screen_pd.rds"))
pd = pd[!is.infinite(pd$pseudotime),]

cell_gene = readRDS(paste0(work_path, "/processing/all_TF_screen_sgRNA/cell_gene.rds"))
df = cell_gene %>% select(cell, gene) %>% unique() %>%
    left_join(pd[,c("cell","pseudotime")], by = "cell") %>% filter(!is.na(pseudotime))
df_num = df %>% group_by(gene) %>% tally() %>% filter(n > 50)
df = df[df$gene %in% as.vector(df_num$gene),]
df$pseudotime_bin = cut(df$pseudotime, breaks = 10, labels = FALSE)

df_x = df %>% group_by(gene, pseudotime_bin) %>% tally() %>%
    dcast(gene ~ pseudotime_bin)
rownames(df_x) = df_x[,1]
df_x = as.matrix(df_x[,-1])
df_x[is.na(df_x)] = 0
df_x = df_x/apply(df_x, 1, sum)

target_list = rownames(df_x)
target_list = target_list[target_list != "NONTARGETING"]

library(philentropy)
res = NULL
for(target_i in target_list){
    print(target_i)
    mat = rbind(df_x[target_i,], rep(0.1, 10))
    res = rbind(res, data.frame(gene = target_i,
                                js_dist = distance(mat, method = "jensen-shannon")))
}
rownames(res) = NULL
res = res[order(res$js_dist, decreasing=T),]

df_x_sub = df_x[res$gene[1:100],]
row_name = NULL
for(i in 1:10){
    for(j in rownames(df_x_sub)){
        if(!j %in% row_name){
            if(which.max(df_x_sub[j,]) == i){
                row_name = c(row_name, j)
            }
        }
    }
}
df_x_sub = df_x[c(row_name),]
write.table(row_name, paste0(work_path, "/analysis/revision/pseudotime/heatmap_big_screen.row_name.txt"), row.names=F, col.names=F, sep="\t", quote=F)

library("gplots")
library(RColorBrewer)
library(viridis)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf(paste0(work_path, "/analysis/revision/pseudotime/heatmap_big_screen.pdf"), 10, 10)
heatmap.2(as.matrix(df_x_sub), 
          col=Colors, 
          #col=viridis(100), 
          scale="row", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(10,5))
dev.off()







