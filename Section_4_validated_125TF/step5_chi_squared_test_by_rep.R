
#####################################################
### Performing chi-squared test (split by replicates)
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tf
### validated_125TF.cell_gene.rds

pd = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned_pd.rds"))
pd$cell = gsub("EB21_CRISPRcut_125TF_GEx_", "", as.vector(pd$cell_id))
cell_gene = readRDS("validated_125TF.cell_gene.rds")


#############
### REP-1 ###
#############

cell_gene = subset(cell_gene, replicate == "rep1")

##########################################################
### Step-1: filtering out targets with low number of cells

### this is the target level analysis, ignoring the heterogenity between guides within the same target

target_cell_num = cell_gene %>% filter(gene != "NONTARGETING") %>% 
    select(cell, gene) %>% unique() %>%
    group_by(gene) %>% tally()

NTC_guide_num = cell_gene %>% filter(gene == "NONTARGETING") %>% 
    select(cell, barcode) %>% unique() %>%
    group_by(barcode) %>% tally()

### filtering out targets with less than 50 cells

target_cell_num = target_cell_num %>% filter(n >= 50)
print(nrow(target_cell_num)) ### n = 113

### filtering out NTCs with less than 20 cells

NTC_guide_num = NTC_guide_num %>% filter(n >= 20)
print(nrow(NTC_guide_num)) ### n = 74


##############################################################
### Step-2: performing chi-squared test on each target vs. ntc

mat = cell_gene %>% group_by(gene, celltype) %>% 
    tally() %>% dcast(gene~celltype, fill = 0)
rownames(mat) = as.vector(mat$gene)
mat = mat[,-1]
mat = as.matrix(mat)

res = NULL
for(i in 1:nrow(target_cell_num)){
    target_i = as.vector(target_cell_num$gene)[i]
    contingency_table = matrix(c(mat[target_i, ], mat["NONTARGETING", ]), nrow = 2, byrow = TRUE)
    chi_sq_test = chisq.test(contingency_table)
    res = c(res, chi_sq_test$p.value)
}
res = data.frame(gene = as.vector(target_cell_num$gene),
                 pval = res)


######################################################################
### Step-3: performing chi-squared test on each two NTC guides vs. ntc

mat = cell_gene %>% filter(gene == "NONTARGETING") %>% group_by(barcode, celltype) %>% 
    tally() %>% dcast(barcode~celltype, fill = 0)
rownames(mat) = as.vector(mat$barcode)
mat = mat[,-1]
mat = as.matrix(mat)

mat_sum = apply(mat, 2, sum)

res_NTC = NULL
for(iter in 1:1){
    print(iter)
    
    each_group = 3
    group_num = floor(nrow(NTC_guide_num)/each_group) + 1
    set.seed(1234)
    NTC_guide_num$group = sample(rep(1:group_num, each = each_group))[1:nrow(NTC_guide_num)]
    
    for(i in 1:group_num){
        barcode_include = as.vector(NTC_guide_num$barcode[NTC_guide_num$group == i])
        
        mat_sub = cell_gene %>% filter(barcode %in% barcode_include) %>% 
            #            slice_sample(n = 100) %>%
            group_by(barcode, celltype) %>% 
            tally() %>% dcast(barcode~celltype, fill = 0)
        rownames(mat_sub) = as.vector(mat_sub$barcode)
        mat_sub = apply(as.matrix(mat_sub[,-1]), 2, sum)
        mat_sub = mat_sub[names(mat_sum)]
        mat_sub[is.na(mat_sub)] = 0
        
        contingency_table = matrix(c(mat_sub, mat_sum), nrow = 2, byrow = TRUE)
        chi_sq_test = chisq.test(contingency_table)
        res_NTC = c(res_NTC, chi_sq_test$p.value)
    }
}

res_sig = res[res$pval < quantile(res_NTC, 0.05),] ### 1.19e-10
res_sig = res_sig[order(res_sig$pval),]
saveRDS(list(res, res_NTC, res_sig), paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/chi_squared_test_rep1.rds"))


############################
### Step-4: making a QQ plot

res_combine = data.frame(pval = c(res$pval, res_NTC),
                         log10pval = c(-log10(res$pval), -log10(res_NTC)),
                         group = c(rep("TF", nrow(res)), rep("NTC", length(res_NTC))))
print(sort(res_combine$log10pval[!is.infinite((res_combine$log10pval))]))
res_combine$log10pval[is.infinite((res_combine$log10pval))] = 165

qlog10 = function(dat){
    # assuming "nn_de" is a dataframe where each row contains one of your p-values
    n <- length(dat)
    
    # Generate quantiles from a normal distribution
    quantiles <- qnorm(seq(0, 1, length.out = n + 1)[-1])
    
    # Rescale quantiles to fit between 0.00001 and 1
    perfect_values <- pmin(pmax(pnorm(quantiles), 1e-5), 1)
    
    # Apply -log10() function to each value
    log_dis <- -log10(perfect_values)
    
    expected_p <- log_dis[order(-log_dis)]
    
    return(expected_p)
}

dat_1 = subset(res_combine, group == "TF")
dat_1 = dat_1[order(dat_1$pval),]
dat_1$expected_p = qlog10(dat_1$pval)

dat_2 = subset(res_combine, group == "NTC")
dat_2 = dat_2[order(dat_2$pval),]
dat_2$expected_p = qlog10(dat_2$pval)

p = ggplot() + 
    geom_point(data = dat_1, aes(x = expected_p, y = log10pval), color = "blue") + 
    geom_point(data = dat_2, aes(x = expected_p, y = log10pval), color = "red") + 
    geom_hline(yintercept = -log10(quantile(res_NTC, 0.05)), linetype="dotted") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/chi_squared_test_QQplot_rep1.pdf"), 4, 4)
print(p)
dev.off()





#############
### REP-2 ###
#############

cell_gene = subset(cell_gene, replicate == "rep2")

##########################################################
### Step-5: filtering out targets with low number of cells

### this is the target level analysis, ignoring the heterogenity between guides within the same target

target_cell_num = cell_gene %>% filter(gene != "NONTARGETING") %>% 
    select(cell, gene) %>% unique() %>%
    group_by(gene) %>% tally()

NTC_guide_num = cell_gene %>% filter(gene == "NONTARGETING") %>% 
    select(cell, barcode) %>% unique() %>%
    group_by(barcode) %>% tally()

### filtering out targets with less than 50 cells

target_cell_num = target_cell_num %>% filter(n >= 50)
print(nrow(target_cell_num)) ### n = 113

### filtering out NTCs with less than 20 cells

NTC_guide_num = NTC_guide_num %>% filter(n >= 20)
print(nrow(NTC_guide_num)) ### n = 76




##############################################################
### Step-6: performing chi-squared test on each target vs. ntc

mat = cell_gene %>% group_by(gene, celltype) %>% 
    tally() %>% dcast(gene~celltype, fill = 0)
rownames(mat) = as.vector(mat$gene)
mat = mat[,-1]
mat = as.matrix(mat)

res = NULL
for(i in 1:nrow(target_cell_num)){
    target_i = as.vector(target_cell_num$gene)[i]
    contingency_table = matrix(c(mat[target_i, ], mat["NONTARGETING", ]), nrow = 2, byrow = TRUE)
    chi_sq_test = chisq.test(contingency_table)
    res = c(res, chi_sq_test$p.value)
}
res = data.frame(gene = as.vector(target_cell_num$gene),
                 pval = res)


######################################################################
### Step-7: performing chi-squared test on each two NTC guides vs. ntc

mat = cell_gene %>% filter(gene == "NONTARGETING") %>% group_by(barcode, celltype) %>% 
    tally() %>% dcast(barcode~celltype, fill = 0)
rownames(mat) = as.vector(mat$barcode)
mat = mat[,-1]
mat = as.matrix(mat)

mat_sum = apply(mat, 2, sum)

res_NTC = NULL
for(iter in 1:1){
    print(iter)
    
    each_group = 3
    group_num = floor(nrow(NTC_guide_num)/each_group) + 1
    set.seed(1234)
    NTC_guide_num$group = sample(rep(1:group_num, each = each_group))[1:nrow(NTC_guide_num)]
    
    for(i in 1:group_num){
        barcode_include = as.vector(NTC_guide_num$barcode[NTC_guide_num$group == i])
        
        mat_sub = cell_gene %>% filter(barcode %in% barcode_include) %>% 
            #            slice_sample(n = 100) %>%
            group_by(barcode, celltype) %>% 
            tally() %>% dcast(barcode~celltype, fill = 0)
        rownames(mat_sub) = as.vector(mat_sub$barcode)
        mat_sub = apply(as.matrix(mat_sub[,-1]), 2, sum)
        mat_sub = mat_sub[names(mat_sum)]
        mat_sub[is.na(mat_sub)] = 0
        
        contingency_table = matrix(c(mat_sub, mat_sum), nrow = 2, byrow = TRUE)
        chi_sq_test = chisq.test(contingency_table)
        res_NTC = c(res_NTC, chi_sq_test$p.value)
    }
}

res_sig = res[res$pval < quantile(res_NTC, 0.05),] ### 1.96e-50
res_sig = res_sig[order(res_sig$pval),]
saveRDS(list(res, res_NTC, res_sig), paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/chi_squared_test_rep2.rds"))


############################
### Step-8: making a QQ plot

res_combine = data.frame(pval = c(res$pval, res_NTC),
                         log10pval = c(-log10(res$pval), -log10(res_NTC)),
                         group = c(rep("TF", nrow(res)), rep("NTC", length(res_NTC))))
print(sort(res_combine$log10pval[!is.infinite((res_combine$log10pval))]))
res_combine$log10pval[is.infinite((res_combine$log10pval))] = 90

qlog10 = function(dat){
    # assuming "nn_de" is a dataframe where each row contains one of your p-values
    n <- length(dat)
    
    # Generate quantiles from a normal distribution
    quantiles <- qnorm(seq(0, 1, length.out = n + 1)[-1])
    
    # Rescale quantiles to fit between 0.00001 and 1
    perfect_values <- pmin(pmax(pnorm(quantiles), 1e-5), 1)
    
    # Apply -log10() function to each value
    log_dis <- -log10(perfect_values)
    
    expected_p <- log_dis[order(-log_dis)]
    
    return(expected_p)
}

dat_1 = subset(res_combine, group == "TF")
dat_1 = dat_1[order(dat_1$pval),]
dat_1$expected_p = qlog10(dat_1$pval)

dat_2 = subset(res_combine, group == "NTC")
dat_2 = dat_2[order(dat_2$pval),]
dat_2$expected_p = qlog10(dat_2$pval)

p = ggplot() + 
    geom_point(data = dat_1, aes(x = expected_p, y = log10pval), color = "blue") + 
    geom_point(data = dat_2, aes(x = expected_p, y = log10pval), color = "red") + 
    geom_hline(yintercept = -log10(quantile(res_NTC, 0.05)), linetype="dotted") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
pdf(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/plot/chi_squared_test_QQplot_rep2.pdf"), 4, 4)
print(p)
dev.off()

