
###########################################################
### Gene expression for those targets which are knocked out
### Chengxiang Qiu
### Feb-20, 2025

source("help_script.R")
mouse_gene <- read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

run_list = c("rep1_EB21_CRISPRi_sorted_gRNA",
             "rep2_EB21_CRISPRcut_sorted_gRNA",
             "rep2_EB21_CRISPRi_sorted_gRNA")

result = list()

for(cnt in 1:3){
    run_id = run_list[cnt]
    print(run_id)
    
    run_id_2 = gsub("_gRNA","",run_id)
    
    cell_barcode = read.table(paste0(work_path, "/processing/pilot_screen_sgRNA/", run_id, "/cell_gene.txt"), header=T, as.is=T, sep="\t")
    cell_barcode$cell = paste0("pilot_screen_GEx_", as.vector(cell_barcode$cell))
    x = as.vector(cell_barcode$gene)
    x[cell_barcode$gene %in% c("NTC.control-random-chr3-116819213-116820756-1",
                               "NTC.control-random-chr3-116819213-116820756-2",
                               "NTC.control-random-chr5-31878254-31878457-3",
                               "NTC.control-random-chr5-31878254-31878457-4",
                               "NTC.control-random-chr7-114416771-114416984-5",
                               "NTC.control-random-chr7-114416771-114416984-6",
                               "NTC.CRISPRCUT-NONTARGETING",
                               "NTC.CRISPRI-NONTARGETING")] = "NTC"
    cell_barcode$gene = as.vector(x)
    
    ### keep cells with single barcodes
    cell_barcode_uniq = cell_barcode %>% group_by(cell) %>% tally() %>% filter(n == 1)
    cell_barcode = cell_barcode[cell_barcode$cell %in% as.vector(cell_barcode_uniq$cell),]
    cell_gene_num = cell_barcode %>% group_by(gene) %>% tally() %>% filter(n >= 20)
    cell_barcode = cell_barcode[cell_barcode$gene %in% as.vector(cell_gene_num$gene),]
    
    gene_list = unique(cell_barcode$gene)
    gene_list = gene_list[gene_list != "NTC"]
    
    mouse_gene_sub = mouse_gene[mouse_gene$gene_short_name %in% gene_list,]
    rownames(mouse_gene_sub) = as.vector(mouse_gene_sub$gene_short_name)
    mouse_gene_sub = mouse_gene_sub[gene_list,]
    
    ### read gene expression data
    obj = readRDS(paste0(work_path, "/processing/pilot_screen_GEx/", run_id_2, "/obj_pilot_screen_GEx_", run_id_2, ".rds"))
    exp = GetAssayData(obj, slot = "counts")
    exp = exp[,as.vector(cell_barcode$cell)]
    exp = t(t(exp) / colSums(exp)) * 10000
    exp@x = log(exp@x + 1)
    
    exp = exp[as.vector(mouse_gene_sub$gene_ID),]
    rownames(exp) = as.vector(mouse_gene_sub$gene_short_name)
    
    exp_aggr = NULL
    for(gene_i in c(gene_list)){
        exp_aggr = cbind(exp_aggr, Matrix::rowMeans(exp[,cell_barcode$gene == gene_i]))
    }
    colnames(exp_aggr) = c(gene_list)
    
    result[[run_id]] = data.frame(gene = gene_list,
                                  gene_exp = as.vector(diag(exp_aggr)),
                                  NTC_exp = Matrix::rowMeans(exp[,cell_barcode$gene == "NTC"]))
    
}

saveRDS(result, paste0(work_path, "/processing/pilot_screen_sgRNA/Compare_TF_exp_to_NTC.rds"))


### locally making scatter plot

result = readRDS("~/work/sam_tf/processing/pilot_screen_sgRNA/Compare_TF_exp_to_NTC.rds")

run_list = c("rep1_EB21_CRISPRi_sorted_gRNA",
             "rep2_EB21_CRISPRcut_sorted_gRNA",
             "rep2_EB21_CRISPRi_sorted_gRNA")

cnt = 3
run_id = run_list[[cnt]]
df = result[[run_id]]

gene_group = rep("developmental TFs", nrow(df))
gene_group[df$gene %in% c("Dnmt1","Eed","Ezh2","Kmt2a","Kmt2d")] = "chromatin modifiers"
gene_group[df$gene %in% c("Ctcf","Gabpa","Nrf1")] = "ES TFs"
gene_group[df$gene %in% c("NTC")] = "NTCs"
df$gene_group = factor(gene_group, levels = c("NTCs","chromatin modifiers","ES TFs","developmental TFs"))

axis_lim = max(max(df$gene_exp), max(df$NTC_exp)) + 0.01

p = ggplot(df, aes(NTC_exp, gene_exp, label = gene)) +
    geom_point() + ggrepel::geom_text_repel(aes(color = gene_group)) +
    xlim(-0.01, axis_lim) +
    ylim(-0.01, axis_lim) +
    geom_abline(intercept = 0) +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = pilot_gene_group_color_code) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 








