
print("Loading Monocle3 and Seurat")
suppressMessages(library(monocle3))
suppressMessages(library(Seurat))

print("Loading packages for regular data analysis, e.g. dplyr")
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))

print("Loading packages for plotting, e.g. ggplot2")
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(gridExtra))
suppressMessages(library(viridis))
suppressMessages(library(gplots))


### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tf
### Please let me know if you cannot find them.
### CX Qiu (cxqiu@uw.edu)

######################################
### color plate used in the manuscript

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

TF_screen_color_code = c("T" = "#939598",
                         "Lhx1" = "#00aeef",
                         "Hand1" = "#00a651",
                         "Gata6" = "#73732a",
                         "Hes5" = "#fbb040",
                         "Six3" = "#f15a29",
                         "Lmo2" = "#81d2e8",
                         "Hhex" = "#f1eb59",
                         "Carm1" = "#ed1c24",
                         "NONTARGETING" = "#662d91")

mEB_time_course_color_code = c("mEB_time_course_GEx_day7"  = "#0d0887",
                               "mEB_time_course_GEx_day14" = "#9c179e",
                               "mEB_time_course_GEx_day21" = "#ed7953")

pilot_screen_GEx_color_code = c("mEB_time_course_GEx_day21" = "#ed7953",
                                "pilot_screen_GEx_rep1_EB21_CRISPRcut_sorted" = "#ac60c3",
                                "pilot_screen_GEx_rep1_EB21_CRISPRcut_unsorted" = "#5fa65c",
                                "pilot_screen_GEx_rep1_EB21_CRISPRi_sorted" = "#ca5c8a",
                                "pilot_screen_GEx_rep1_EB21_CRISPRi_unsorted" = "#6a86ce",
                                "pilot_screen_GEx_rep2_EB21_CRISPRcut_sorted" = "#c14845",
                                "pilot_screen_GEx_rep2_EB21_CRISPRi_sorted" = "#ae943c")

pilot_gene_group_color_code = c("chromatin modifiers" = "#a46cb7",
                                "ES TFs" = "#4866d3",
                                "developmental TFs" = "#5dc598",
                                "NTCs" = "#808080")

###########################################
### published datasets used for integration

jax_celltype_color_code = c("Epithelium" = "#9c476b",
                            "Mesoderm" = "#5dbb4d",
                            "Endothelium" = "#7459c6",
                            "Cardiomyocytes" = "#a8b439",
                            "Megakaryocytes" = "#c845a1",
                            "CNS_neurons" = "#46873f",
                            "Hepatocytes" = "#bc76d5",
                            "Olfactory_sensory_neurons" = "#d59a35",
                            "Ependymal_cells" = "#6c80c9",
                            "Lung_and_airway" = "#d05130",
                            "Neuroectoderm_and_glia" = "#4eacd7",
                            "Intestine" = "#73732a",
                            "Neural_crest_PNS_glia" = "#5dc598",
                            "Primitive_erythroid" = "#d681b2",
                            "Neural_crest_PNS_neurons" = "#9cb26a",
                            "Muscle_cells" = "#d3756f",
                            "White_blood_cells" = "#379884",
                            "Eye_and_other" = "#9f5f2b",
                            "Definitive_erythroid" = "#d69f6a")

pijuan_celltype_color_code = c("Epiblast" = "#5dbb4d",
                               "ExE ectoderm" = "#5dc54f",
                               "ExVE" = "#d59a35",
                               "ExE endoderm" = "#93b730",
                               "Nascent mesoderm" = "#d05130",
                               "PGC" = "#9c476b",
                               "Rostral neurectoderm" = "#7459c6",
                               "Def. endoderm" = "#bc76d5",
                               "Mesenchyme" = "#4eacd7",
                               "Blood progenitors 2" = "#85bf71",
                               "Gut" = "#4866d3",
                               "Erythroid3" = "#cabd39",
                               "Notochord" = "#46873f",
                               "Caudal Mesoderm" =              "#dc9535",
                               "Paraxial mesoderm" =            "#446bae",
                               "Somitic mesoderm" =             "#d4552a",
                               "Allantois" =                    "#4cacdb",
                               "Forebrain/Midbrain/Hindbrain" = "#d4424c",
                               "Cardiomyocytes" =               "#73732a",
                               "Neural crest" =                 "#b4397f",
                               "Primitive Streak" =             "#347e47",
                               "EmVE" =                         "#da84cf",
                               "Parietal endoderm" = "#6c80c9",
                               "Visceral endoderm" = "#b282d5",
                               "Anterior Primitive Streak" = "#aa9732",
                               "Mixed mesoderm" = "#894e98",
                               "Surface ectoderm" = "#c845a1",
                               "Haematoendothelial progenitors" = "#88a0e5",
                               "Blood progenitors 1" = "#d681b2",
                               "ExE mesoderm" = "#67b88c",
                               "Caudal epiblast" = "#459939",
                               "Erythroid2" = "#2b7f63",
                               "Intermediate mesoderm" = "#dd7f9e",
                               "Pharyngeal mesoderm" = "#547739",
                               "Caudal neurectoderm" = "#ea8e72",
                               "Erythroid1" = "#9cb26a",
                               "Endothelium" = "#5dc598",
                               "Spinal cord" = "#bead6c",
                               "NMP" = "#b98554")

Liberali_celltype_color_code = c("Naive pluripotency" = "#9c476b",
                                 "Exiting naive pluripotency" = "#9cb26a",
                                 "Zscan4+ Artefact" = "#7459c6",
                                 "Epiblast" = "#a8b439",
                                 "Pre-somitic mesoderm" = "#c845a1",
                                 "Somite differentiation front" = "#46873f",
                                 "Ectopic pluripotency" = "#9f5f2b",
                                 "Caudal mesoderm" = "#d59a35",
                                 "Gut" = "#6c80c9",
                                 "Paraxial mesoderm" = "#d05130",
                                 "Hemogenic endothelium" = "#5dc598",
                                 "Somite" = "#73732a",
                                 "Neuromesodermal progenitors" = "#4eacd7",
                                 "Cd63+ ectoderm-like artefact" = "#d681b2",
                                 "Epiblast/primitive streak" = "#5dbb4d",
                                 "Primitive streak" = "#d3756f",
                                 "Caudal epiblast/primitive streak" = "#379884",
                                 "Anterior primitive streak/Def. endoderm" = "#bc76d5",
                                 "Caudal epiblast" = "#d69f6a")

TLS_celltype_color_code = c("PCGLC" = "#9c476b",
                            "NeuralTube2" = "#9cb26a",
                            "NeuralTube1" = "#7459c6",
                            "NMPs" = "#a8b439",
                            "pPSM" = "#c845a1",
                            "aPSM" = "#46873f",
                            "Somite-1" = "#9f5f2b",
                            "Somite0" = "#d59a35",
                            "Somite" = "#6c80c9",
                            "SomiteDermo" = "#d05130",
                            "SomiteSclero" = "#5dc598",
                            "Endothelial" = "#73732a",
                            "Endoderm" = "#4eacd7",
                            "Unknown" = "#d681b2")

Rosen_celltype_color_code = c("Primitive Streak"          = "#dc3c6e",
                              "Early Nascent Mesoderm"    = "#5cc151",
                              "Anterior Primitive Streak" = "#b15ecf",
                              "Early Posterior PSM"       = "#9bb837",
                              "Mature Endoderm"           = "#6a61cd",
                              "Caudal Epiblast"           = "#c9a93a",
                              "PGCs"                      = "#6576c1",
                              "Early Neurectoderm"        = "#df852f",
                              "Epiblast"                  = "#5fa1d8",
                              "Caudal Neurectoderm"       = "#d54637",
                              "Late Posterior PSM"        = "#42c3c2",
                              "Early Spinal Cord"         = "#d149a2",
                              "Late Neurectoderm"         = "#4f923a",
                              "Notochord"                 = "#d28fd0",
                              "Late Nascent Mesoderm"     = "#66c388",
                              "Cardiopharyngeal Mesoderm" = "#974f8a",
                              "NMPs"                      = "#747c29",
                              "Mature Somites"            = "#df7a89",
                              "Somites"                   = "#3f926d",
                              "Early Anterior PSM"        = "#a44657",
                              "Head Mesoderm"             = "#476b31",
                              "Late Spinal Cord"          = "#a8532d",
                              "Late Anterior PSM"         = "#aeb06a",
                              "Endothelium"               = "#e1956d",
                              "Neurons"                   = "#957030")

Van_celltype_color_code = c("Cardiac"               = "#be73a7",
                            "Paraxial MD"           = "#5eb956",
                            "Differentiated somite" = "#c152b8",
                            "Somite"                = "#b0b33e",
                            "Differentiation front" = "#7b66cd",
                            "PSM"                   = "#d88c38",
                            "NMPs"                  = "#678ccc",
                            "Spinal cord"           = "#ce4c33",
                            "Mesenchyme"            = "#50b9a3",
                            "Endothelium"           = "#cc4270",
                            "Allantois"             = "#558140",
                            "PGC-like/ExE EcD"      = "#c56f61",
                            "Endoderm"              = "#9c8240")


#####################################################
### Function: doing regular analysis using Seurat ###
#####################################################

doClusterSeurat <- function(obj, 
                            nfeatures = 2500, 
                            resolution = 1, 
                            k.filter = 200, 
                            doClustering = TRUE, 
                            n_dim = 30, 
                            min.dist = 0.3, 
                            n.components = 2){
    
    if(length(table(obj$group))!=1){
        
        obj.list <- SplitObject(object = obj, split.by = "group")
        
        for (i in 1:length(x = obj.list)) {
            obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
            obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], 
                                                  selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
        }
        
        reference.list <- obj.list[names(table(obj$group))]
        obj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:n_dim, k.filter = k.filter)
        
        obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:n_dim)
        
        DefaultAssay(object = obj.integrated) <- "integrated"
        
        obj <- obj.integrated 
        
    } else {
        
        obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
        
    }
    
    obj <- ScaleData(object = obj, verbose = FALSE)
    obj <- RunPCA(object = obj, npcs = n_dim, verbose = FALSE)
    
    if (doClustering == TRUE){
        obj <- FindNeighbors(object = obj, dims = 1:n_dim, reduction = "pca")
        obj <- FindClusters(object = obj, resolution = resolution)
    }
    
    obj <- RunUMAP(object = obj, reduction = "pca", dims = 1:n_dim, min.dist = min.dist, n.components = n.components)
    
    return(obj)
}

