
#################################################
### Processing sgRNA data in the pilot experiment
### Chengxiang Qiu
### Feb-20, 2025

###############################
### Step-1: processing 10x data

### rep1 ###

# EB21_CRISPRcut_unsorted_gRNA 
# EB21_CRISPRcut_sorted_1_gRNA 
# EB21_CRISPRcut_sorted_2_gRNA
# EB21_CRISPRi_unsorted_gRNA 
# EB21_CRISPRi_sorted_1_gRNA 
# EB21_CRISPRi_sorted_2_gRNA 


module load cellranger/7.2.0

cellranger count \
--fastqs /net/shendure/vol8/projects/ajh24/ajh24/proj/2018eb_factor_screen/data/reads/2018_10_22_eb_pilot_screen/cellranger_outputs.guides/HWGC5BGX5/outs/fastq_path \
--nosecondary \
--localcores 8 \
--include-introns true \
--chemistry SC3Pv2 \
--sample EB21_CRISPRcut_unsorted_gRNA \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_sgRNA/rep1_EB21_CRISPRcut_unsorted_gRNA \
--id EB21_CRISPRcut_unsorted_gRNA \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0

cellranger count \
--fastqs /net/shendure/vol8/projects/ajh24/ajh24/proj/2018eb_factor_screen/data/reads/2018_10_22_eb_pilot_screen/cellranger_outputs.guides/HWGC5BGX5/outs/fastq_path \
--nosecondary \
--localcores 8 \
--include-introns true \
--chemistry SC3Pv2 \
--sample EB21_CRISPRi_unsorted_gRNA \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_sgRNA/rep1_EB21_CRISPRi_unsorted_gRNA \
--id EB21_CRISPRi_unsorted_gRNA \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0

cellranger count \
--fastqs /net/shendure/vol8/projects/ajh24/ajh24/proj/2018eb_factor_screen/data/reads/2018_10_22_eb_pilot_screen/cellranger_outputs.guides/HWGC5BGX5/outs/fastq_path \
--nosecondary \
--localcores 8 \
--include-introns true \
--chemistry SC3Pv2 \
--sample EB21_CRISPRcut_sorted_1_gRNA,EB21_CRISPRcut_sorted_2_gRNA \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_sgRNA/rep1_EB21_CRISPRcut_sorted_gRNA \
--id EB21_CRISPRcut_sorted_gRNA \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0

cellranger count \
--fastqs /net/shendure/vol8/projects/ajh24/ajh24/proj/2018eb_factor_screen/data/reads/2018_10_22_eb_pilot_screen/cellranger_outputs.guides/HWGC5BGX5/outs/fastq_path \
--nosecondary \
--localcores 8 \
--include-introns true \
--chemistry SC3Pv2 \
--sample EB21_CRISPRi_sorted_1_gRNA,EB21_CRISPRi_sorted_2_gRNA \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_sgRNA/rep1_EB21_CRISPRi_sorted_gRNA \
--id EB21_CRISPRi_sorted_gRNA \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0


### rep2 ###

# EB21_CRISPRcut_sorted_1_gRNA 
# EB21_CRISPRcut_sorted_2_gRNA 
# EB21_CRISPRcut_sorted_3_gRNA 
# EB21_CRISPRi_sorted_1_gRNA 
# EB21_CRISPRi_sorted_2_gRNA 
# EB21_CRISPRi_sorted_3_gRNA


module load cellranger/7.2.0

cellranger count \
--fastqs /net/shendure/vol10/projects/Samuel/nobackup/10X/2019_EB21_pilot_screen/HYJYBGXB_gRNA_1/outs/fastq_path/,/net/shendure/vol10/projects/Samuel/nobackup/10X/2019_EB21_pilot_screen/HJCFKBGXB_gRNA_2/outs/fastq_path/ \
--nosecondary \
--localcores 8 \
--include-introns true \
--chemistry SC3Pv3 \
--sample EB21_CRISPRcut_sorted_1_gRNA,EB21_CRISPRcut_sorted_1_gRNA_2,EB21_CRISPRcut_sorted_2_gRNA,EB21_CRISPRcut_sorted_2_gRNA_2,EB21_CRISPRcut_sorted_3_gRNA,EB21_CRISPRcut_sorted_3_gRNA_2 \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_sgRNA/rep2_EB21_CRISPRcut_sorted_gRNA \
--id EB21_CRISPRcut_sorted_gRNA \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0

cellranger count \
--fastqs /net/shendure/vol10/projects/Samuel/nobackup/10X/2019_EB21_pilot_screen/HYJYBGXB_gRNA_1/outs/fastq_path/,/net/shendure/vol10/projects/Samuel/nobackup/10X/2019_EB21_pilot_screen/HJCFKBGXB_gRNA_2/outs/fastq_path/ \
--nosecondary \
--localcores 8 \
--include-introns true \
--chemistry SC3Pv3 \
--sample EB21_CRISPRi_sorted_1_gRNA,EB21_CRISPRi_sorted_1_gRNA_2,EB21_CRISPRi_sorted_2_gRNA,EB21_CRISPRi_sorted_2_gRNA_2,EB21_CRISPRi_sorted_3_gRNA,EB21_CRISPRi_sorted_3_gRNA_2 \
--output-dir /net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_sgRNA/rep2_EB21_CRISPRi_sorted_gRNA \
--id EB21_CRISPRi_sorted_gRNA \
--transcriptome /net/shendure/vol10/projects/Samuel/nobackup/10X/refdata-cellranger-mm10-3.0.0


#########################################################################
### Step-2: running get_barcode.py to generate table with cell x barcodes

conda deactivate
module load python/2.7.13
module load pysam/0.15.2

data_path=/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/data/pilot_screen_sgRNA
output_path=/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf/processing/pilot_screen_sgRNA

python /net/gs/vol1/home/cxqiu/work/scripts/sam_tf/get_barcodes_silvia.py \
--input_bams "$data_path"/rep1_EB21_CRISPRcut_unsorted_gRNA/outs/possorted_genome_bam.bam \
-o "$output_path"/rep1_EB21_CRISPRcut_unsorted_gRNA/get_bc_output.txt \
--whitelist "$output_path"/crisprcut.whitelist.txt \
--search_seq GTGGAAAGGACGAAACACCG \
--no_swalign

python /net/gs/vol1/home/cxqiu/work/scripts/sam_tf/get_barcodes_silvia.py \
--input_bams "$data_path"/rep1_EB21_CRISPRcut_sorted_gRNA/outs/possorted_genome_bam.bam \
-o "$output_path"/rep1_EB21_CRISPRcut_sorted_gRNA/get_bc_output.txt \
--whitelist "$output_path"/crisprcut.whitelist.txt \
--search_seq GTGGAAAGGACGAAACACCG \
--no_swalign

python /net/gs/vol1/home/cxqiu/work/scripts/sam_tf/get_barcodes_silvia.py \
--input_bams "$data_path"/rep1_EB21_CRISPRi_unsorted_gRNA/outs/possorted_genome_bam.bam \
-o "$output_path"/rep1_EB21_CRISPRi_unsorted_gRNA/get_bc_output.txt \
--whitelist "$output_path"/crispri.whitelist.txt \
--search_seq GTGGAAAGGACGAAACACCG \
--no_swalign

python /net/gs/vol1/home/cxqiu/work/scripts/sam_tf/get_barcodes_silvia.py \
--input_bams "$data_path"/rep1_EB21_CRISPRi_sorted_gRNA/outs/possorted_genome_bam.bam \
-o "$output_path"/rep1_EB21_CRISPRi_sorted_gRNA/get_bc_output.txt \
--whitelist "$output_path"/crispri.whitelist.txt \
--search_seq GTGGAAAGGACGAAACACCG \
--no_swalign

python /net/gs/vol1/home/cxqiu/work/scripts/sam_tf/get_barcodes_silvia.py \
--input_bams "$data_path"/rep2_EB21_CRISPRcut_sorted_gRNA/outs/possorted_genome_bam.bam \
-o "$output_path"/rep2_EB21_CRISPRcut_sorted_gRNA/get_bc_output.txt \
--whitelist "$output_path"/crisprcut.whitelist.txt \
--search_seq GTGGAAAGGACGAAACACCG \
--no_swalign

python /net/gs/vol1/home/cxqiu/work/scripts/sam_tf/get_barcodes_silvia.py \
--input_bams "$data_path"/rep2_EB21_CRISPRi_sorted_gRNA/outs/possorted_genome_bam.bam \
-o "$output_path"/rep2_EB21_CRISPRi_sorted_gRNA/get_bc_output.txt \
--whitelist "$output_path"/crispri.whitelist.txt \
--search_seq GTGGAAAGGACGAAACACCG \
--no_swalign


#############################################################
### Step-3: running preprocess_cfg_features.R to filter sgRNA

### After processing the GEx, we generate the cell list which was used in the following analysis (as white list)

run_id_list = c("rep1_EB21_CRISPRcut_unsorted",
                "rep1_EB21_CRISPRcut_sorted",
                "rep2_EB21_CRISPRcut_sorted",
                "rep1_EB21_CRISPRi_unsorted",
                "rep1_EB21_CRISPRi_sorted",
                "rep2_EB21_CRISPRi_sorted")
library(Seurat)
work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf"
for(run_id in run_id_list){
    print(run_id)
    obj = readRDS(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/obj_", "pilot_screen_GEx_", run_id, ".rds"))
    x = gsub("pilot_screen_GEx_", "", colnames(obj))
    write.table(x, paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/cell_list_by_manual_filtering.txt"), row.names=F, col.names=F, sep="\t", quote=F)
}


####################################################################
### Step-4: filtering sgRNAs and connecting to TFs based on barcodes

run_id_list = c("rep1_EB21_CRISPRcut_unsorted",
                "rep1_EB21_CRISPRcut_sorted",
                "rep2_EB21_CRISPRcut_sorted",
                "rep1_EB21_CRISPRi_unsorted",
                "rep1_EB21_CRISPRi_sorted",
                "rep2_EB21_CRISPRi_sorted")

for(cnt in 1:length(run_id_list)){
    
    run_id = run_id_list[cnt]
    print(run_id)
    
    cell_white_list = read.table(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/cell_list_by_manual_filtering.txt"), as.is=T, sep="\t")
    colnames(cell_white_list) = "cell"
    
    cell_barcode = read.table(paste0(work_path, "/processing/pilot_screen_sgRNA/", run_id, "_gRNA/get_bc_output.txt"), header=T, as.is=T, sep="\t")
    cell_barcode = cell_barcode[!str_detect(cell_barcode$barcode, 'unprocessed'), ]
    cell_barcode$log2_read_count = log2(cell_barcode$read_count + 1)
    cell_barcode$log2_umi_count = log2(cell_barcode$umi_count + 1)
    
    cell_barcode_sub = subset(cell_barcode, read_count >= 10 & umi_count >= 2)
    p1 = cell_barcode_sub %>% 
        ggplot(aes(log2_read_count, log2_umi_count)) + geom_point()
    p2 = cell_barcode_sub %>%
        ggplot(aes(log2_read_count)) + geom_histogram(bins = 50) 
    p3 = cell_barcode_sub %>%
        ggplot(aes(log2_umi_count)) + geom_histogram(bins = 50) 
    
    pdf(paste0(work_path, "/processing/pilot_screen_sgRNA/", run_id, "_gRNA/bc_umi_read_dis.pdf"), 10, 4)
    grid.arrange(p1, p2, p3, nrow=1, ncol=3) 
    dev.off()
    
}

### After manually reviewing the distribution plots, we set those cutoffs:

log2_umi_cutoff = 2
log2_read_cutoff = 5

df = NULL

for(cnt in 1:length(run_id_list)){
    
    run_id = run_id_list[cnt]
    
    cell_white_list = read.table(paste0(work_path, "/processing/pilot_screen_GEx/", run_id, "/cell_list_by_manual_filtering.txt"), as.is=T, sep="\t")
    colnames(cell_white_list) = "cell"
    
    cell_barcode = read.table(paste0(work_path, "/processing/pilot_screen_sgRNA/", run_id, "_gRNA/get_bc_output.txt"), header=T, as.is=T, sep="\t")
    cell_barcode = cell_barcode[!str_detect(cell_barcode$barcode, 'unprocessed'), ]
    cell_barcode$cell =  paste0(run_id, "_", unlist(lapply(as.vector(cell_barcode$cell), function(x) strsplit(x,"[-]")[[1]][1])))
    cell_barcode$log2_read_count = log2(cell_barcode$read_count + 1)
    cell_barcode$log2_umi_count = log2(cell_barcode$umi_count + 1)
    
    cell_barcode_sub = subset(cell_barcode, log2_umi_count >= log2_umi_cutoff & log2_read_count >= log2_read_cutoff)
    cell_barcode_sub$log2_read_count = cell_barcode_sub$log2_umi_count = NULL
    
    if (str_detect(run_id, 'cut')){
        tf_barcode = read.table(paste0(work_path, "/processing/pilot_screen_sgRNA/crisprcut.plus_ntcs.guide_metadata.txt"), as.is=T)
    } else {
        tf_barcode = read.table(paste0(work_path, "/processing/pilot_screen_sgRNA/crispri.plus_ntcs.guide_metadata.txt"), as.is=T)
    }
    colnames(tf_barcode) = c("gene", "barcode")
    if(sum(! cell_barcode_sub$barcode %in% tf_barcode$barcode) != 0){
        print("some barcode doesn't map to any TF")
    }
    
    cell_barcode_sub = cell_barcode_sub %>% left_join(tf_barcode, by = "barcode")
    cell_barcode_sub = subset(cell_barcode_sub, cell %in% cell_white_list$cell)
    
    #write.table(cell_barcode_sub, paste0(work_path, "/processing/pilot_screen_sgRNA/", run_id, "_gRNA/cell_gene.txt"), row.names=F, quote=F, sep="\t")
    
    a = round(100*length(unique(cell_barcode_sub$cell))/length(unique(cell_white_list$cell)), 2)
    b = cell_barcode_sub %>% group_by(cell) %>% tally() %>% filter(n == 1)
    c = round(100*nrow(b)/length(unique(cell_barcode_sub$cell)), 2)
    print(paste0(run_id, " : ", a, "%", " cells with sgRNA, and of those ", c, "% had only one sgRNA"))
    
    df = rbind(df, data.frame(experiment_id = run_id,
                              total_cell = length(unique(cell_white_list$cell)),
                              barcode_cell = length(unique(cell_barcode_sub$cell)),
                              uniq_barcode_cell = nrow(b)))
}

# "rep1_EB21_CRISPRcut_unsorted : 0.22% cells with sgRNA, and of those 100% had only one sgRNA"
# "rep1_EB21_CRISPRcut_sorted : 27.76% cells with sgRNA, and of those 98.25% had only one sgRNA"
# "rep2_EB21_CRISPRcut_sorted : 81.28% cells with sgRNA, and of those 77.29% had only one sgRNA"
# "rep1_EB21_CRISPRi_unsorted : 0.78% cells with sgRNA, and of those 100% had only one sgRNA"
# "rep1_EB21_CRISPRi_sorted : 43.98% cells with sgRNA, and of those 92.73% had only one sgRNA"
# "rep2_EB21_CRISPRi_sorted : 81.42% cells with sgRNA, and of those 60.47% had only one sgRNA"

df = df[df$experiment_id != "rep1_EB21_CRISPRcut_sorted",]
df$barcode_pct = 100*df$barcode_cell/df$total_cell
df$uniq_barcode_pct = 100*df$uniq_barcode_cell/df$total_cell

df2 = data.frame(experiment_id = c(as.vector(df$experiment_id), as.vector(df$experiment_id)),
                 pct = c(df$barcode_pct, df$uniq_barcode_pct),
                 group = c(rep("barcode", nrow(df)), rep("uniq_barcode", nrow(df))))
df2$experiment_id = factor(df2$experiment_id, levels = rev(c("rep1_EB21_CRISPRcut_unsorted",
                                                           "rep1_EB21_CRISPRi_unsorted",
                                                           "rep2_EB21_CRISPRcut_sorted",
                                                           "rep1_EB21_CRISPRi_sorted",
                                                           "rep2_EB21_CRISPRi_sorted")))
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


###############################
### Step-5: making summary plot

### a pie chart of the plasmid library breakdown (~ 8% ES-essential TFs (n=3), 
### 13% chromatin modifiers (n=5), 60% developmental TFs (n=23), 20% NTCs (n = x))

work_path = "/net/shendure/vol10/projects/cxqiu/nobackup/work/sam_tf"
dat = read.table(paste0(work_path, "/processing/pilot_screen_sgRNA/crisprcut.plus_ntcs.guide_metadata.txt"), as.is=T, sep="\t")
names(dat) = c("gene", "barcode")

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

gene_group = rep("developmental TFs", nrow(dat))
gene_group[dat$gene %in% c("Dnmt1","Eed","Ezh2","Kmt2a","Kmt2d")] = "chromatin modifiers"
gene_group[dat$gene %in% c("Ctcf","Gabpa","Nrf1")] = "ES TFs"
gene_group[dat$gene %in% c("NTC")] = "NTCs"
dat$gene_group = factor(gene_group, levels = c("NTCs","chromatin modifiers","ES TFs","developmental TFs"))

data = dat %>% group_by(gene_group) %>% tally() %>% rename(count = n)
data$gene_group = factor(data$gene_group, levels = c("NTCs","chromatin modifiers","ES TFs","developmental TFs"))

# Compute percentages
data$fraction <- data$count / sum(data$count)

# Compute the cumulative percentages (top of each rectangle)
data$ymax <- cumsum(data$fraction)

# Compute the bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))

# Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

# Compute a good label
data$label <- paste0(data$gene_group, "\n", data$count)

# Make the plot
p = ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=gene_group)) +
    geom_rect() +
    geom_text( x=2, aes(y=labelPosition, label=label, color=gene_group), size=6) + # x here controls label position (inner / outer)
    scale_fill_manual(values = pilot_gene_group_color_code) +
    scale_color_manual(values = pilot_gene_group_color_code) +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) +
    theme_void() +
    theme(legend.position = "none")

pdf("~/share/TF_bc_pie.pdf")
print(p)
dev.off()








