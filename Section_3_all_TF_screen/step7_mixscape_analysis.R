
############################################
### Performing mixscape analysis on mesoderm
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]


obj = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned.rds"))
pd = readRDS(paste0(work_path, "/analysis/all_TF_screen_GEx/obj_aligned_pd.rds"))
pd$cell = gsub("all_TF_screen_GEx_", "", as.vector(pd$cell_id))

cell_gene = readRDS(paste0(work_path, "/processing/all_TF_screen_sgRNA/cell_gene.rds"))
cell_gene = unique(cell_gene[,c("cell", "gene")])
cell_gene_sub = subset(cell_gene, cell %in% as.vector(pd$cell[pd$celltype == "Lateral plate mesoderm"]))
cell_gene_num = cell_gene_sub %>% group_by(cell) %>% tally() %>% filter(n == 1)
cell_gene_sub = subset(cell_gene_sub, cell %in% as.vector(cell_gene_num$cell))
gene_num = table(cell_gene_sub$gene)
cell_gene_sub = subset(cell_gene_sub, gene %in% names(gene_num[gene_num > 50]))
### 3696 cells, across 13 targets and NTC

pd_sub = pd[pd$cell %in% as.vector(cell_gene_sub$cell),]
count_sub = GetAssayData(obj, slot = "count")[,pd$cell %in% as.vector(cell_gene_sub$cell)]
pd_sub_x = pd_sub %>% left_join(cell_gene_sub %>% select(cell, gene), by = "cell")
pd_sub$gene = as.vector(pd_sub_x$gene)

meso = CreateSeuratObject(count_sub, meta.data = pd_sub)
meso = NormalizeData(meso, normalization.method = "LogNormalize", scale.factor = 10000)
meso <- FindVariableFeatures(meso, selection.method = "vst", nfeatures = 2500)
meso <- ScaleData(object = meso, verbose = FALSE)
meso <- RunPCA(object = meso, npcs = 30, verbose = FALSE)

# Setup custom theme for plotting.
custom_theme <- theme(plot.title = element_text(size = 16, hjust = 0.5), legend.key.size = unit(0.7, 
                                                                                                "cm"), legend.text = element_text(size = 14))

# Calculate local perturbation signature.
meso <- CalcPerturbSig(object = meso, assay = "RNA", slot = "data", gd.class = "gene", nt.cell.class = "NONTARGETING", 
                       reduction = "pca", ndims = 30, num.neighbors = 20, new.assay.name = "PRTB")


# Prepare PRTB assay for dimensionality reduction: Normalize data, find variable features and
# center data.
DefaultAssay(object = meso) <- "PRTB"

# Use variable features from RNA assay.
VariableFeatures(object = meso) <- VariableFeatures(object = meso[["RNA"]])
meso <- ScaleData(meso, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
meso <- RunPCA(object = meso, reduction.key = "prtbpca", reduction.name = "prtbpca")

# Run UMAP to visualize clustering in 2-D.
meso <- RunUMAP(object = meso, dims = 1:40, reduction = "prtbpca", reduction.key = "prtbumap", 
                reduction.name = "prtbumap")


# Run mixscape to classify cells based on their perturbation status.
meso <- RunMixscape(object = meso, assay = "PRTB", slot = "scale.data", labels = "gene", nt.class.name = "NONTARGETING", 
                    min.de.genes = 5, iter.num = 10, de.assay = "RNA", verbose = F)


# Remove non-perturbed cells and run LDA to reduce the dimensionality of the data.
Idents(meso) <- "mixscape_class.global"
sub <- subset(meso, idents = c("KO", "NONTARGETING"))

sub <- MixscapeLDA(object = sub, assay = "RNA", pc.assay = "PRTB", labels = "gene", nt.label = "NONTARGETING", 
                   npcs = 10, logfc.threshold = 0.25, verbose = F)

# Use LDA results to run UMAP and visualize cells on 2-D.
sub <- RunUMAP(sub, dims = 1:11, reduction = "lda", reduction.key = "ldaumap", reduction.name = "ldaumap")

# Visualize UMAP clustering results.
Idents(sub) <- "mixscape_class"
sub$mixscape_class <- as.factor(sub$mixscape_class)
p <- DimPlot(sub, reduction = "ldaumap", label = T, repel = T, label.size = 5)

library(scales)
col = setNames(object = hue_pal()(14), nm = levels(sub$mixscape_class))
names(col) <- c(names(col)[1:10], "NONTARGETING", names(col)[12:14])
col[11] <- "grey39"

pdf(paste0(work_path, '/analysis/all_TF_screen_GEx/plot/mixscape_mesoderm_umap.pdf'), width=10, height=5)
p + scale_color_manual(values = col, drop = FALSE) + ylab("UMAP 2") + xlab("UMAP 1") + custom_theme
dev.off()

saveRDS(list(meso, sub), paste0(work_path, '/analysis/all_TF_screen_GEx/mixscape_mesoderm.rds'))

#Knockout Efficiencies
table(sub$mixscape_class)/table(meso$gene)

Asz1 KO      Cbfb KO      Cbx8 KO     Cecr6 KO    Ctnnb1 KO      Gsc2 KO
0.7924528    0.7142857    0.9375000    0.8985507    0.9833333    0.8229167
Hnrnpr KO    Hoxa11 KO      Lhx8 KO  Map3k7cl KO NONTARGETING     Nr4a2 KO
0.7500000    0.8139535    1.0000000    0.9021739    1.0000000    0.9531250
Rfx2 KO     Trp53 KO
0.6986301    0.6911076


### identifying markers for individual target

perturb=unique(meso$mixscape_class[grep('KO',meso$mixscape_class)])

markers = NULL
for (i in 1:length(perturb)){
    print(i)
    per <- perturb[i]
    markers_i <- FindMarkers(sub, ident.1 = per, ident.2='NONTARGETING', min.pct = 0.25)
    markers_i$gene = gsub(" KO", "", per)
    markers_i$gene_ID = rownames(markers_i)
    markers = rbind(markers, markers_i)
}
rownames(markers) = NULL
markers = markers %>% left_join(mouse_gene[,c("gene_ID","gene_short_name")], by = "gene_ID")
saveRDS(markers, paste0(work_path, '/analysis/all_TF_screen_GEx/mixscape_mesoderm_markers.rds'))

markers_filter = markers %>% filter(avg_log2FC > 0.25, p_val_adj < 0.001)
markers_filter = markers_filter[,c("gene", "avg_log2FC", "p_val_adj", "gene_short_name")]
write.table(markers_filter, "~/share/mixscape_mesoderm_markers_sig.txt", row.names=F, sep="\t", quote=F)


### perform GO analysis for individual target
library(topGO)
target_list = unique(markers$gene)

res = NULL
for(i in target_list){
    print(i)
    markers_i = subset(markers, gene == i)
    cluster = markers_i %>% filter(avg_log2FC > 0.25, p_val_adj < 0.001)
    
    geneList <- ifelse(unique(as.vector(markers_i$gene_short_name)) %in% as.vector(cluster$gene_short_name), 1, 0)
    names(geneList) <- unique(as.vector(markers_i$gene_short_name)) #BETTER
    
    # Create topGOdata object
    GOdata <- new("topGOdata",
                  ontology = "BP", 
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
    resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
    res_i = GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
    res_i$gene = i
    
    res = rbind(res, res_i)
}
saveRDS(res, paste0(work_path, '/analysis/all_TF_screen_GEx/mixscape_mesoderm_markers_GO.rds'))

i = "Trp53"
res_sub = subset(res, gene == i)
res_sub$log10_pval = -log10(as.numeric(res_sub$Fisher))
res_sub$Term = factor(res_sub$Term, levels = rev(as.vector(res_sub$Term)))

p<-ggplot(data=res_sub, aes(x=Term, y=log10_pval)) +
    geom_bar(stat="identity") + 
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    coord_flip() + 
    ggsave(paste0(work_path, '/analysis/all_TF_screen_GEx/plot/mixscape_mesoderm_markers_GO_', i ,'.pdf'))






