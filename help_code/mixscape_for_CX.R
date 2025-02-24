#R4.0.0
module load pcre2/10.35
module load R/4.0.0

#For mixscape
library(remotes)
remotes::install_github('satijalab/seurat', ref = 'mixscape')

library(Seurat)
#library(SeuratData) #not available
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
options(stringsAsFactors=FALSE)

##load EB seurat objects and metadata (with gRNA info)
seurat <- readRDS('/net/shendure/vol10/projects/silvia/EBs/scaled/nobackup/cut/all_cells_transferred_celltype.RDS')
meta <- read.table('/net/shendure/vol10/projects/silvia/EBs/scaled/nobackup/cut/reseq_monocle_metadata.txt', sep='\t', header=TRUE)

seurat <- UpdateSeuratObject(seurat)

seurat@meta.data = meta[rownames(seurat@meta.data),]

#select cells with 1 guide only
meso <- subset(seurat, cells=rownames(seurat@meta.data[seurat$guide_count==1,]))

#select 1 celltype and only guides with at least 50 cells, here: mesoderm
meso <- subset(meso, cells=rownames(meso@meta.data[meso$celltype=='pharyngeal mesoderm',]))

freq_guides <- names(table(meso$gene))[table(meso$gene)>50] #14 + ntc
meso <- subset(meso, cells=rownames(meso@meta.data[meso$gene %in% freq_guides,])) #3364 cells


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

col = setNames(object = hue_pal()(15), nm = levels(sub$mixscape_class))
names(col) <- c(names(col)[1:11], "NONTARGETING", names(col)[13:15])
col[12] <- "grey39"

pdf('/net/shendure/vol10/projects/silvia/EBs/scaled/nobackup/cut/mixscape_mesoderm_umap.pdf', width=10, height=5)
p + scale_color_manual(values = col, drop = FALSE) + ylab("UMAP 2") + xlab("UMAP 1") + custom_theme
dev.off()




p2 <- DimPlot(meso, reduction = "umap", label = F, repel = T, label.size = 5)
col2 = setNames(object = hue_pal()(15), nm = levels(meso$gene))
names(col2) <- c( "NONTARGETING", names(col2)[2:15])
col2[1] <- "grey39"

col2 <- col
names(col2) <- str_replace(names(col2),' KO','')

pdf('/net/shendure/vol10/projects/silvia/EBs/scaled/nobackup/cut/mesoderm_umap_no_mixscape_no_label.pdf', width=10, height=5)
p2 + scale_color_manual(values = col2, drop = FALSE) + ylab("UMAP 2") + xlab("UMAP 1") + custom_theme
dev.off()


#Knockout Efficiencies
table(sub$mixscape_class)/table(meso$gene)

    Agap3 KO      Asz1 KO      Cbfb KO      Cbx8 KO     Cecr6 KO    Ctnnb1 KO 
   0.7500000    0.8854167    0.8983051    0.9642857    0.8571429    0.8571429 
     Gsc2 KO    Hnrnpr KO    Hoxa11 KO      Lhx8 KO  Map3k7cl KO NONTARGETING 
   0.6153846    0.8072289    0.8658537    1.0000000    0.9010989    1.0000000 
    Nr4a2 KO      Rfx2 KO     Trp53 KO 
   0.9649123    0.7384615    0.7567084 


lhx8.markers <- FindMarkers(sub, ident.1 = 'Lhx8 KO', ident.2='NONTARGETING', min.pct = 0.25)

cbx8.markers <- FindMarkers(sub, ident.1 = 'Cbx8 KO', ident.2='NONTARGETING', min.pct = 0.25)

rfx2.markers <- FindMarkers(sub, ident.1 = 'Rfx2 KO', ident.2='NONTARGETING', min.pct = 0.25)

trp53.markers <- FindMarkers(sub, ident.1 = 'Trp53 KO', ident.2='NONTARGETING', min.pct = 0.25)

ctnnb.markers <- FindMarkers(sub, ident.1 = 'Ctnnb1 KO', ident.2='NONTARGETING', min.pct = 0.25)
hoxa11.markers <- FindMarkers(sub, ident.1 = 'Hoxa11 KO', ident.2='NONTARGETING', min.pct = 0.25)
nr4a2.markers <- FindMarkers(sub, ident.1 = 'Nr4a2 KO', ident.2='NONTARGETING', min.pct = 0.25)

perturb=unique(meso$mixscape_class[grep('KO',meso$mixscape_class)])

for (i in 1:length(perturb)){
	per <- perturb[i]
	markers <- FindMarkers(sub, ident.1 = per, ident.2='NONTARGETING', min.pct = 0.25)
	write.table(markers, paste0('/net/shendure/vol10/projects/silvia/EBs/scaled/nobackup/cut/mixscape_mesoderm_markers_',per,'.txt'),quote=F, row.names=T, sep='\t')
	}

cluster <- lhx8.markers[lhx8.markers$avg_logFC > 0.25 & lhx8.markers$p_val_adj < 0.001, ]
cluster <- cbx8.markers[cbx8.markers$avg_logFC > 0.25 & cbx8.markers$p_val_adj < 0.001, ]
cluster <- trp53.markers[trp53.markers$avg_logFC > 0.25 & trp53.markers$p_val_adj < 0.001, ]
cluster <- rfx2.markers[rfx2.markers$avg_logFC > 0.25 & rfx2.markers$p_val_adj < 0.001, ]
cluster <- nr4a2.markers[nr4a2.markers$avg_logFC > 0.25 & nr4a2.markers$p_val_adj < 0.001, ]
cluster <- hoxa11.markers[hoxa11.markers$avg_logFC > 0.25 & hoxa11.markers$p_val_adj < 0.001, ]
cluster <- ctnnb.markers[ctnnb.markers$avg_logFC > 0.25 & ctnnb.markers$p_val_adj < 0.001, ]

expr <- rownames(cluster)
# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(unique(rownames(lhx8.markers)) %in% expr, 1, 0)
names(geneList) <- unique(rownames(lhx8.markers)) #BETTER


# Create topGOdata object
    GOdata <- new("topGOdata",
        ontology = "BP", 
        allGenes = geneList,
        geneSelectionFun = function(x)(x == 1), annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
    resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)

#lhx8 bone development (heart, kidney, lung)
#trp53 dna damage, antiviral
#cbx8 forelimb, mesenchyme
#rfx2 metabolism

sum(unique(sub$gene) %in% list_filt[,1])




#Try with select NTC as conrol


meso.ntc <- subset(meso, cells=rownames(seurat@meta.data[seurat$gene=='NONTARGETING',]))




#Try random NTCs 

set.seed(1234)
Ctrl_1 <- sample(rownames(meso.ntc@meta.data), 91)
set.seed(1234)
Ctrl_2 <- sample(rownames(meso.ntc@meta.data)[!rownames(meso.ntc@meta.data) %in% Ctrl_1], 66)

meso@meta.data[rownames(meso@meta.data) %in% Ctrl_1,13] <- 'Ctrl_1'
meso@meta.data[rownames(meso@meta.data) %in% Ctrl_2,13] <- 'Ctrl_2'

#run same code as above, then

# Visualize UMAP clustering results.
Idents(sub) <- "mixscape_class"
sub$mixscape_class <- as.factor(sub$mixscape_class)
p <- DimPlot(sub, reduction = "ldaumap", label = T, repel = T, label.size = 5)

col = setNames(object = hue_pal()(15), nm = levels(sub$mixscape_class))
names(col) <- c(names(col)[1:11], "NONTARGETING", names(col)[13:15])
col[12] <- "grey39"

pdf('/net/shendure/vol10/projects/silvia/EBs/scaled/nobackup/cut/mixscape_mesoderm_umap_random_ntcs.pdf', width=10, height=5)
p + scale_color_manual(values = col, drop = FALSE) + ylab("UMAP 2") + xlab("UMAP 1") + custom_theme
dev.off()
#100% are NP -> this works
