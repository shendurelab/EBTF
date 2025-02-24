
###########################################################
### Processing sgRNA data in the validated 125TF experiment
### Chengxiang Qiu
### Feb-20, 2025

###############################################
### Step-1, processing sgRNA using 10x pipeline

#!/bin/bash
INPUTFILES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32)
INPUTFILENAME="${INPUTFILES[$SGE_TASK_ID - 1]}"
echo $INPUTFILENAME

module load cellranger/7.2.0

cellranger count \
--fastqs /net/shendure/vol9/nobackup/sam_tf/SD_SR_125TF_EB_HTO/SD_SR_125TF_EB_HTO_done/,/net/shendure/vol10/projects/silvia/EBs/scaled_125TFs/nobackup/cellranger_gRNA/AAAG727HV_gRNA/outs/fastq_path/AAAG727HV/EB_gRNA_"$INPUTFILENAME"/ \
--nosecondary \
--localcores 8 \
--include-introns true \
--chemistry SC3Pv3 \
--sample EB_gRNA_"$INPUTFILENAME" \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/EB21_CRISPRcut_125TF_sgRNA/EB_"$INPUTFILENAME"_gRNA \
--id EB_"$INPUTFILENAME"_gRNA \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0


### running get_barcode.py to generate table with cell x barcodes ###

INPUTFILES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32)
INPUTFILENAME="${INPUTFILES[$SGE_TASK_ID - 1]}"
echo $INPUTFILENAME

conda deactivate
module load python/2.7.13
module load pysam/0.15.2

data_path=/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/EB21_CRISPRcut_125TF_sgRNA
output_path=/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/processing/EB21_CRISPRcut_125TF_sgRNA

python /net/gs/vol1/home/cxqiu/work/scripts/sam_tf/get_barcodes_silvia.py \
--input_bams "$data_path"/EB_"$INPUTFILENAME"_gRNA/outs/possorted_genome_bam.bam \
-o "$output_path"/EB_"$INPUTFILENAME"_gRNA/get_bc_output.txt \
--whitelist "$output_path"/whitelist_cut_125TF_gRNA_get_barcode.txt \
--search_seq GTGGAAAGGACGAAACACCG \
--no_swalign



#############################################################
### Step-2: running preprocess_cfg_features.R to filter sgRNA

### After processing the GEx, we generate the cell list which was used in the following analysis (as white list)

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf"
pd = readRDS(paste0(work_path, "/analysis/EB21_CRISPRcut_125TF_GEx/obj_aligned_pd.rds"))
x = data.frame(cell_barcode = unlist(lapply(as.vector(pd$cell_id), function(x) strsplit(x,"[_]")[[1]][7])),
               cell_group = gsub("EB21_CRISPRcut_125TF_GEx_", "", pd$experiment_id))
write.table(x, paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/cell_list_by_manual_filtering.txt"), row.names=F, col.names=F, sep="\t", quote=F)


####################################################################
### Step-3: filtering sgRNAs and connecting to TFs based on barcodes

run_id_list = paste0("EB_", c(1:32))

getCutoff = function(x){
    y = 1.25*x - 5
    return(y)
}

log2_umi_cutoff = 2

cell_white_list = read.table(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/cell_list_by_manual_filtering.txt"), as.is=T, sep="\t")
colnames(cell_white_list) = c("cell","experiment_id")
cell_white_list$cell = paste0(cell_white_list$experiment_id, "_", cell_white_list$cell)

tf_barcode = read.table(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/guide_crispr_cut_scaled_followup_125TF_metadata.txt"), as.is=T)
colnames(tf_barcode) = c("gene", "barcode")

tmp = as.vector(tf_barcode$gene)
tmp[tf_barcode$gene %in% c("control_random_chr3:116819213-116820756_1",
                           "control_random_chr3:116819213-116820756_2",
                           "control_random_chr5:31878254-31878457_3",
                           "control_random_chr5:31878254-31878457_4",
                           "control_random_chr7:114416771-114416984_5",
                           "control_random_chr7:114416771-114416984_6")] = "NONTARGETING"
tf_barcode$gene = as.vector(tmp)

df = NULL

for(cnt in 1:length(run_id_list)){
    
    run_id = run_id_list[cnt]
    print(run_id)
    
    cell_barcode = read.table(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/", run_id, "_gRNA/get_bc_output.txt"), header=T, as.is=T, sep="\t")
    cell_barcode = cell_barcode[!str_detect(cell_barcode$barcode, 'unprocessed'), ]
    cell_barcode$cell =  paste0(run_id, "_", unlist(lapply(as.vector(cell_barcode$cell), function(x) strsplit(x,"[-]")[[1]][1])))
    cell_barcode$log2_read_count = log2(cell_barcode$read_count + 1)
    cell_barcode$log2_umi_count = log2(cell_barcode$umi_count + 1)
    cell_barcode$id = paste0(cell_barcode$cell, "_", cell_barcode$barcode)
    
    ### before filtering low-quality gRNAs
    cell_barcode_sub = subset(cell_barcode, cell %in% cell_white_list$cell)
    a = round(100*length(unique(cell_barcode_sub$cell))/length(unique(cell_white_list$cell[cell_white_list$experiment_id == run_id])), 2)
    b = cell_barcode_sub %>% group_by(cell) %>% tally() %>% filter(n == 1)
    c = round(100*nrow(b)/length(unique(cell_barcode_sub$cell)), 2)

    df = rbind(df, data.frame(experiment_id = run_id,
                              total_cell = length(unique(cell_white_list$cell[cell_white_list$experiment_id == run_id])),
                              barcode_cell = length(unique(cell_barcode_sub$cell)),
                              uniq_barcode_cell = nrow(b)))
    
    ### after filtering low-quality gRNAs by UMI count and read count
    cell_barcode$tmp = getCutoff(cell_barcode$log2_read_count)
    cell_barcode_x = subset(cell_barcode, log2_umi_count < tmp)
    cell_barcode_x = subset(cell_barcode_x, log2_umi_count >= log2_umi_cutoff)
    
    cell_barcode_sub = subset(cell_barcode, read_count >= 10 & umi_count >= 2)
    cell_barcode_sub$filter = if_else(cell_barcode_sub$id %in% as.vector(cell_barcode_x$id), "yes", "no")
    
    p1 = cell_barcode_sub %>% 
        ggplot(aes(log2_read_count, log2_umi_count, color = filter)) + geom_point() +
        scale_color_manual(values=c("yes" = "red", "no" = "grey")) + theme(legend.position="none") +
        geom_abline(intercept = -5, slope = 1.25, color = "black") + geom_hline(yintercept = log2_umi_cutoff, color = "black")
    p2 = cell_barcode_sub %>%
        ggplot(aes(log2_read_count)) + geom_histogram(bins = 50) 
    p3 = cell_barcode_sub %>%
        ggplot(aes(log2_umi_count)) + geom_histogram(bins = 50) 
    
    pdf(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/", run_id, "_gRNA/bc_umi_read_dis.pdf"), 10, 4)
    grid.arrange(p1, p2, p3, nrow=1, ncol=3) 
    dev.off()
    
}

df$barcode_pct = 100*df$barcode_cell/df$total_cell
df$uniq_barcode_pct = 100*df$uniq_barcode_cell/df$total_cell

saveRDS(df, paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/pct_barcode_per_lanes_before_filtering.rds"))

print(sum(df$barcode_cell)/sum(df$total_cell))
print(sum(df$uniq_barcode_cell)/sum(df$barcode_cell))
### Before filtering out low-quality sgRNAs
### Across the 32 samples, 87.4% of cells with at least one sgRNA, of which 23.6% has only one single sgRNA

### After manually reviewing the distribution plots, we decide those cutoffs:

getCutoff = function(x){
    y = 1.25*x - 5
    return(y)
}

log2_umi_cutoff = 2

df = NULL

for(cnt in 1:length(run_id_list)){
    
    run_id = run_id_list[cnt]
    
    cell_barcode = read.table(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/", run_id, "_gRNA/get_bc_output.txt"), header=T, as.is=T, sep="\t")
    cell_barcode = cell_barcode[!str_detect(cell_barcode$barcode, 'unprocessed'), ]
    cell_barcode$cell =  paste0(run_id, "_", unlist(lapply(as.vector(cell_barcode$cell), function(x) strsplit(x,"[-]")[[1]][1])))
    cell_barcode$log2_read_count = log2(cell_barcode$read_count + 1)
    cell_barcode$log2_umi_count = log2(cell_barcode$umi_count + 1)
    
    if(sum(! cell_barcode_sub$barcode %in% tf_barcode$barcode) != 0){
        print("some barcode doesn't map to any TF")
    }
    
    cell_barcode$tmp = getCutoff(cell_barcode$log2_read_count)
    cell_barcode_sub = subset(cell_barcode, log2_umi_count < tmp)
    cell_barcode_sub = subset(cell_barcode_sub, log2_umi_count >= log2_umi_cutoff)
    
    cell_barcode_sub$log2_read_count = cell_barcode_sub$log2_umi_count = cell_barcode_sub$tmp = NULL
    
    cell_barcode_sub = cell_barcode_sub %>% left_join(tf_barcode, by = "barcode")
    cell_barcode_sub = subset(cell_barcode_sub, cell %in% cell_white_list$cell)
    
    write.table(cell_barcode_sub, paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/", run_id, "_gRNA/cell_gene.txt"), row.names=F, quote=F, sep="\t")
    
    a = round(100*length(unique(cell_barcode_sub$cell))/length(unique(cell_white_list$cell[cell_white_list$experiment_id == run_id])), 2)
    b = cell_barcode_sub %>% group_by(cell) %>% tally() %>% filter(n == 1)
    c = round(100*nrow(b)/length(unique(cell_barcode_sub$cell)), 2)
    print(paste0(run_id, " : ", a, "%", " cells with at least one sgRNA, and of those ", c, "% had only one sgRNA"))
    
    df = rbind(df, data.frame(experiment_id = run_id,
                              total_cell = length(unique(cell_white_list$cell[cell_white_list$experiment_id == run_id])),
                              barcode_cell = length(unique(cell_barcode_sub$cell)),
                              uniq_barcode_cell = nrow(b)))
}

[1] "EB_1 : 56.45% cells with at least one sgRNA, and of those 95.06% had only one sgRNA"
[1] "EB_2 : 57.85% cells with at least one sgRNA, and of those 95.51% had only one sgRNA"
[1] "EB_3 : 64.77% cells with at least one sgRNA, and of those 90.94% had only one sgRNA"
[1] "EB_4 : 63.62% cells with at least one sgRNA, and of those 91.93% had only one sgRNA"
[1] "EB_5 : 62.6% cells with at least one sgRNA, and of those 90.69% had only one sgRNA"
[1] "EB_6 : 66.14% cells with at least one sgRNA, and of those 89.8% had only one sgRNA"
[1] "EB_7 : 59.68% cells with at least one sgRNA, and of those 92.51% had only one sgRNA"
[1] "EB_8 : 68.78% cells with at least one sgRNA, and of those 90.46% had only one sgRNA"
[1] "EB_9 : 46.47% cells with at least one sgRNA, and of those 83.29% had only one sgRNA"
[1] "EB_10 : 63.3% cells with at least one sgRNA, and of those 64.92% had only one sgRNA"
[1] "EB_11 : 60.93% cells with at least one sgRNA, and of those 68.57% had only one sgRNA"
[1] "EB_12 : 55.85% cells with at least one sgRNA, and of those 72.99% had only one sgRNA"
[1] "EB_13 : 58% cells with at least one sgRNA, and of those 71.8% had only one sgRNA"
[1] "EB_14 : 51.51% cells with at least one sgRNA, and of those 77.98% had only one sgRNA"
[1] "EB_15 : 65.5% cells with at least one sgRNA, and of those 65.32% had only one sgRNA"
[1] "EB_16 : 47.02% cells with at least one sgRNA, and of those 80.68% had only one sgRNA"
[1] "EB_17 : 80.56% cells with at least one sgRNA, and of those 78.77% had only one sgRNA"
[1] "EB_18 : 79.93% cells with at least one sgRNA, and of those 80.71% had only one sgRNA"
[1] "EB_19 : 79.63% cells with at least one sgRNA, and of those 80.94% had only one sgRNA"
[1] "EB_20 : 80.76% cells with at least one sgRNA, and of those 80.12% had only one sgRNA"
[1] "EB_21 : 80.16% cells with at least one sgRNA, and of those 80.55% had only one sgRNA"
[1] "EB_22 : 79.94% cells with at least one sgRNA, and of those 81.28% had only one sgRNA"
[1] "EB_23 : 80.1% cells with at least one sgRNA, and of those 80.9% had only one sgRNA"
[1] "EB_24 : 81.24% cells with at least one sgRNA, and of those 82.36% had only one sgRNA"
[1] "EB_25 : 77.37% cells with at least one sgRNA, and of those 83.47% had only one sgRNA"
[1] "EB_26 : 79.64% cells with at least one sgRNA, and of those 81.11% had only one sgRNA"
[1] "EB_27 : 80.31% cells with at least one sgRNA, and of those 79.27% had only one sgRNA"
[1] "EB_28 : 81.6% cells with at least one sgRNA, and of those 81.39% had only one sgRNA"
[1] "EB_29 : 80.75% cells with at least one sgRNA, and of those 81.33% had only one sgRNA"
[1] "EB_30 : 79.71% cells with at least one sgRNA, and of those 80.4% had only one sgRNA"
[1] "EB_31 : 78.37% cells with at least one sgRNA, and of those 80.42% had only one sgRNA"
[1] "EB_32 : 80.43% cells with at least one sgRNA, and of those 80.16% had only one sgRNA"


df$barcode_pct = 100*df$barcode_cell/df$total_cell
df$uniq_barcode_pct = 100*df$uniq_barcode_cell/df$total_cell

saveRDS(df, paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/pct_barcode_per_lanes.rds"))

print(sum(df$barcode_cell)/sum(df$total_cell))
print(sum(df$uniq_barcode_cell)/sum(df$barcode_cell))
### After filtering out low-quality sgRNAs
### Across the 32 samples, 68.5% of cells with at least one sgRNA, of which 79.0% has only one single sgRNA

df2 = data.frame(experiment_id = c(as.vector(df$experiment_id), as.vector(df$experiment_id)),
                 pct = c(df$barcode_pct, df$uniq_barcode_pct),
                 group = c(rep("barcode", nrow(df)), rep("uniq_barcode", nrow(df))))
df2$experiment_id = factor(df2$experiment_id, levels = rev(paste0("EB_", c(1:32))))
df2$group = factor(df2$group, levels = c("uniq_barcode", "barcode"))

p = ggplot(data=df2, aes(x=experiment_id, y=pct, fill = group)) +
    geom_bar(stat="identity", position=position_dodge()) + coord_flip() +
    labs(x = "Experiment", y = "%") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set2") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

pdf(paste0("~/share/barcode_pct.pdf"), 4, 5)
print(p)
dev.off()


########################################################
### Step-4: create a big file for the following analysis 

run_id_list = paste0("EB_", c(1:32))

dat = NULL
for(cnt in 1:length(run_id_list)){
    print(cnt); run_id = run_id_list[cnt]
    dat_i = read.table(paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/", run_id, "_gRNA/cell_gene.txt"),header=T,as.is=T)
    dat = rbind(dat, dat_i)
}

### how to assign TF to each single cell
### to get as many cells assigned for each TF, we don't just include those cells with single gRNA assigned
### i.e. a cell with Gata6_NTC would be counted as Gata6 only, whereas Gata6_Ctcf would be counted as both Gata6 and Ctcf

dat_uniq = dat %>% select(cell, gene) %>% unique() %>%
    group_by(cell) %>% tally() %>% filter(n == 1)
dat_1 = dat[dat$cell %in% as.vector(dat_uniq$cell),]
dat_2 = dat[!dat$cell %in% as.vector(dat_uniq$cell),]
dat_2 = subset(dat_2, gene != "NONTARGETING")

dat_x = rbind(dat_1, dat_2)

saveRDS(dat_x, paste0(work_path, "/processing/EB21_CRISPRcut_125TF_sgRNA/cell_gene.rds"))

### how many cells per target
dat_x$experiment_id = paste0("EB_", unlist(lapply(as.vector(dat_x$cell), function(x) strsplit(x,"[_]")[[1]][2])))
dat_x$replicate = rep("rep1", nrow(dat_x))
dat_x$replicate[dat_x$experiment_id %in% paste0("EB_", c(17:32))] = "rep2"


x = dat_x %>% select(cell, gene, replicate) %>% unique() %>% filter(replicate == "rep1") %>% 
    filter(gene != "NONTARGETING") %>% group_by(gene) %>% tally()
summary(x$n)
### median 342, mean 557, across 125 targets

x = dat_x %>% select(cell, gene, replicate) %>% unique() %>% filter(replicate == "rep2") %>% 
    filter(gene != "NONTARGETING") %>% group_by(gene) %>% tally()
summary(x$n)
### median 413, mean 707, across 125 targets



