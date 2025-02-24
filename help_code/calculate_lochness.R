target_lochNESS <- function(obj, tf) {
    safe_column_scale = function(bmat, scale_vector) {
        bmat@x <- bmat@x / rep.int(scale_vector, diff(bmat@p))
        return(bmat)
    }
    safe_column_multiply = function(bmat, scale_vector) {
        bmat@x <- bmat@x * rep.int(scale_vector, diff(bmat@p))
        return(bmat)
    }
    kadj = round(0.5 * sqrt(ncol(obj)))
    # lochness
    # adding a "genotype" column that is "KO" for cbcs assigned a Cdx2 gRNA,
    # and "WT" for all other cells
    if (tf == "NTC") {
        obj@meta.data$genotype <- ifelse(obj@meta.data$target == "NTC", "KO", "WT")
    } else {
        obj@meta.data$genotype <- ifelse(grepl(paste0('(^|_)', tf, '(_|$)'), obj@meta.data$target) == TRUE, "KO", "WT")
    }
    
    # MT scores
    mt_counts = as.data.frame(table(obj$genotype))
    mutant_mask = ifelse(obj@meta.data$genotype=="WT", 0, 1)
    mutant_neighbors = Matrix::rowSums(safe_column_multiply(obj@graphs$tf.all.peaks_nn, mutant_mask))
    pca_mutant_score = (mutant_neighbors / kadj) / (1 - mt_counts$Freq[mt_counts$Var1=='WT'] / length(mutant_neighbors)) - 1
    
    # plotting
    mutant_score = pca_mutant_score
    df_umap = as.data.frame(obj[['TFallumap']]@cell.embeddings)
    
    df_umap$raw_mt_score = mutant_score
    # df_umap$mt_score = scale(mutant_score, center = TRUE, scale = FALSE)
    
    bg1 = df_umap[obj$genotype=="WT",]
    dt1 = df_umap[obj$genotype!="WT",]
    
    df_umap$genotype <- obj@meta.data$genotype[match(rownames(df_umap), rownames(obj@meta.data))]
    df_umap
}


