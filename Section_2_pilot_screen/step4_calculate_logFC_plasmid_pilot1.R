

#################################################################
### Calculate logFC between sgRNA from plasmid and mEBs (Pilot-1)
### Chengxiang Qiu
### Feb-20, 2025


source("help_script.R")
mouse_gene <- read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

work_path = "Your_work_path"

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/sam_tf

plasmid_cut_file = "CRISPRcut_pilot_1_sgRNAs_from_plasmid_ESC.txt"
plasmid_i_file = "CRISPRi_pilot_1_sgRNAs_from_plasmid_ESC_EB_ind.txt"


#####################################
### Step-1: read plasmid and ESC file

plasmid_cut = read.table(plasmid_cut_file, as.is=T)
plasmid_cut$gene = unlist(lapply(as.vector(plasmid_cut$sgRNAs), function(x) strsplit(x,"[_]")[[1]][1]))
plasmid_cut_1 = plasmid_cut[plasmid_cut$sgRNAs %in% c("NTC_i_chr3:116819213-116820756_1",
                                                      "NTC_i_chr3:116819213-116820756_2",
                                                      "NTC_i_chr5:31878254-31878457_3",
                                                      "NTC_i_chr5:31878254-31878457_4",
                                                      "NTC_i_chr7:114416771-114416984_5",
                                                      "NTC_i_chr7:114416771-114416984_6"),]
plasmid_cut_1$barcode = gsub("NTC_i_","",as.vector(plasmid_cut_1$sgRNAs))
plasmid_cut_2 = plasmid_cut[!plasmid_cut$sgRNAs %in% c("NTC_i_chr3:116819213-116820756_1",
                                                      "NTC_i_chr3:116819213-116820756_2",
                                                      "NTC_i_chr5:31878254-31878457_3",
                                                      "NTC_i_chr5:31878254-31878457_4",
                                                      "NTC_i_chr7:114416771-114416984_5",
                                                      "NTC_i_chr7:114416771-114416984_6"),]
plasmid_cut_2$barcode = unlist(lapply(as.vector(plasmid_cut_2$sgRNAs), function(x) strsplit(x,"[.]")[[1]][2])) 
plasmid_cut_2$barcode = unlist(lapply(as.vector(plasmid_cut_2$barcode), function(x) substr(x, 2, nchar(x) - 1))) 
plasmid_cut = rbind(plasmid_cut_1, plasmid_cut_2)
plasmid_cut$prop_plasmid = plasmid_cut$plasmid_cut/sum(plasmid_cut$plasmid_cut)
plasmid_cut$prop_ESC = plasmid_cut$ES_tr_cut/sum(plasmid_cut$ES_tr_cut)
plasmid_cut$barcode = gsub("[:|_]", "-", plasmid_cut$barcode)
plasmid_cut$barcode = paste0(plasmid_cut$gene, ":", plasmid_cut$barcode)

plasmid_i = read.table(plasmid_i_file, as.is=T)
plasmid_i$plasmid_i = plasmid_i$plasmid_1_i + plasmid_i$plasmid_2_i
plasmid_i$gene = unlist(lapply(as.vector(plasmid_i$sgRNAs), function(x) strsplit(x,"[_]")[[1]][1]))
plasmid_i_1 = plasmid_i[plasmid_i$sgRNAs %in% c("NTC_i_chr3:116819213-116820756_1",
                                                "NTC_i_chr3:116819213-116820756_2",
                                                "NTC_i_chr5:31878254-31878457_3",
                                                "NTC_i_chr5:31878254-31878457_4",
                                                "NTC_i_chr7:114416771-114416984_5",
                                                "NTC_i_chr7:114416771-114416984_6"),]
plasmid_i_1$barcode = gsub("NTC_i_","",as.vector(plasmid_i_1$sgRNAs))
plasmid_i_2 = plasmid_i[!plasmid_i$sgRNAs %in% c("NTC_i_chr3:116819213-116820756_1",
                                                 "NTC_i_chr3:116819213-116820756_2",
                                                 "NTC_i_chr5:31878254-31878457_3",
                                                 "NTC_i_chr5:31878254-31878457_4",
                                                 "NTC_i_chr7:114416771-114416984_5",
                                                 "NTC_i_chr7:114416771-114416984_6"),]
plasmid_i_2$barcode = unlist(lapply(as.vector(plasmid_i_2$sgRNAs), function(x) strsplit(x,"[.]")[[1]][2])) 
plasmid_i_2$barcode = unlist(lapply(as.vector(plasmid_i_2$barcode), function(x) substr(x, 2, nchar(x)))) 
plasmid_i = rbind(plasmid_i_1, plasmid_i_2)
plasmid_i$prop_plasmid = plasmid_i$plasmid_i/sum(plasmid_i$plasmid_i)
plasmid_i$prop_ESC = plasmid_i$EB_ind_i/sum(plasmid_i$EB_ind_i)
plasmid_i$barcode = gsub("[:|_]", "-", plasmid_i$barcode)
plasmid_i$barcode = paste0(plasmid_i$gene, ":", plasmid_i$barcode)


##############################################
### Step-2: calculate the FC and make the plot

run_list = c("rep1_EB21_CRISPRcut_sorted_gRNA", 
             "rep1_EB21_CRISPRi_sorted_gRNA")

res = list()
res_orig = list()

for(run_id in run_list){
    
    print(run_id)
    
    dat = read.table(paste0(work_path, "/processing/pilot_screen_sgRNA/", run_id, "/cell_gene.txt"), header=T, as.is=T, sep="\t")
    
    dat_1 = dat[dat$gene %in% c("NTC.control-random-chr3-116819213-116820756-1",
                                "NTC.control-random-chr3-116819213-116820756-2",
                                "NTC.control-random-chr5-31878254-31878457-3",
                                "NTC.control-random-chr5-31878254-31878457-4",
                                "NTC.control-random-chr7-114416771-114416984-5",
                                "NTC.control-random-chr7-114416771-114416984-6"),]
    dat_1$barcode = gsub("NTC.control-random-", "", as.vector(dat_1$gene))
    dat_2 = dat[!dat$gene %in% c("NTC.control-random-chr3-116819213-116820756-1",
                                "NTC.control-random-chr3-116819213-116820756-2",
                                "NTC.control-random-chr5-31878254-31878457-3",
                                "NTC.control-random-chr5-31878254-31878457-4",
                                "NTC.control-random-chr7-114416771-114416984-5",
                                "NTC.control-random-chr7-114416771-114416984-6"),]
    dat = rbind(dat_1, dat_2)
    
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
    dat$barcode = paste0(dat$gene, ":", dat$barcode)
    
    if(run_id %in% c("rep1_EB21_CRISPRcut_sorted_gRNA","rep2_EB21_CRISPRcut_sorted_gRNA")){
        if(sum(!dat$barcode %in% plasmid_cut$barcode) != 0){
            print("WRONG")
        }
    } else {
        if(sum(!dat$barcode %in% plasmid_i$barcode) != 0){
            print("WRONG")
        }
    }
    
    ### i.e. a cell with Gata6_NTC would be counted as Gata6 only, whereas Gata6_Ctcf would be counted as both Gata6 and Ctcf
    
    dat_uniq = dat %>% group_by(cell) %>% tally() %>% filter(n == 1)
    dat_1 = dat[dat$cell %in% as.vector(dat_uniq$cell),]
    
    #dat_mul = dat %>% group_by(cell) %>% tally() %>% filter(n > 1)
    #dat_2 = dat[dat$cell %in% as.vector(dat_mul$cell) & dat$gene != "NTC",]
    
    #dat = rbind(dat_1, dat_2)
    
    ### I found that using cells with single sgRNA assigned looks better, so giving up with the above approach (maybe better for the later big screen data)
    dat = dat_1
    
    if(run_id %in% c("rep1_EB21_CRISPRcut_sorted_gRNA","rep2_EB21_CRISPRcut_sorted_gRNA")){
        df1 = dat %>% group_by(barcode) %>% tally() %>% rename(count_EB = n) %>% 
            mutate(prop_EB = count_EB/sum(count_EB)) %>% left_join(plasmid_cut, by = "barcode") %>% 
            mutate(log2_fc = log2(prop_EB/prop_plasmid), comparing = "EB") %>% as.data.frame()
        df2 = plasmid_cut %>% filter(barcode %in% df1$barcode) %>%
            mutate(log2_fc = log2(prop_ESC/prop_plasmid), comparing = "ESC") %>% as.data.frame()
        df = rbind(df1[,c("gene", "log2_fc", "comparing")],
                   df2[,c("gene", "log2_fc", "comparing")])
    } else {
        df1 = dat %>% group_by(barcode) %>% tally() %>% rename(count_EB = n) %>% 
            mutate(prop_EB = count_EB/sum(count_EB)) %>% left_join(plasmid_i, by = "barcode") %>% 
            mutate(log2_fc = log2(prop_EB/prop_plasmid), comparing = "EB") %>% as.data.frame()
        df2 = plasmid_i %>% filter(barcode %in% df1$barcode) %>%
            mutate(log2_fc = log2(prop_ESC/prop_plasmid), comparing = "ESC") %>% as.data.frame()
        df = rbind(df1[,c("gene", "log2_fc", "comparing")],
                   df2[,c("gene", "log2_fc", "comparing")])
    }
    
    df = df[!is.infinite(df$log2_fc),]
    
    res_orig[[run_id]] = df
    
    df = df %>% group_by(gene, comparing) %>% 
        summarize(mean_log2_fc = mean(log2_fc), sd_log2_fc = sd(log2_fc), count = n()) %>% as.data.frame()
    df$se_log2_fc = df$sd_log2_fc/sqrt(df$count)
    
    gene_group = rep("developmental TFs", nrow(df))
    gene_group[df$gene %in% c("Dnmt1","Eed","Ezh2","Kmt2a","Kmt2d")] = "chromatin modifiers"
    gene_group[df$gene %in% c("Ctcf","Gabpa","Nrf1")] = "ES TFs"
    gene_group[df$gene %in% c("NTC")] = "NTCs"
    df$gene_group = factor(gene_group, levels = c("NTCs","chromatin modifiers","ES TFs","developmental TFs"))
    
    df$id = paste0(df$comparing, ": ", df$gene)
    df$comparing = factor(df$comparing, levels = c("ESC", "EB"))
    df = df[order(df$comparing, df$gene_group, df$gene),]
    df$id = factor(df$id, levels = as.vector(df$id))
    
    res[[run_id]] = df
    
    p = ggplot(data=df, aes(x=id, y=mean_log2_fc, fill = gene_group)) +
        geom_bar(stat="identity") +
        geom_errorbar(aes(ymin = mean_log2_fc - sd_log2_fc, ymax = mean_log2_fc + sd_log2_fc), width = 0.4) +
        labs(x = "", y = "Log2 gRNA FC to plasmid library", title = run_id) +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = pilot_gene_group_color_code) +
        theme(axis.text.x = element_text(color="black", angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(color="black")) +
        geom_vline(xintercept = 32.5, linetype="dashed", color = "grey50")
    
    pdf(paste0(work_path, "/processing/pilot_screen_sgRNA/", run_id, "/LogFC_EB_vs_plasmid.pdf"), 8, 4)
    print(p)
    dev.off()
}

saveRDS(res, paste0(work_path, "/processing/pilot_screen_sgRNA/rep1_LogFC_EB_vs_plasmid.rds"))
saveRDS(res_orig, paste0(work_path, "/processing/pilot_screen_sgRNA/rep1_LogFC_EB_vs_plasmid_orig.rds"))


########################################################
### Step-3: showing a summary boxplot across four groups

res_orig = readRDS(paste0(work_path, "/processing/pilot_screen_sgRNA/rep1_LogFC_EB_vs_plasmid_orig.rds"))

run_list = c("rep1_EB21_CRISPRcut_sorted_gRNA", 
             "rep1_EB21_CRISPRi_sorted_gRNA")

for(run_id in run_list){
    print(run_id)
    df = res_orig[[run_id]]
    
    gene_group = rep("developmental TFs", nrow(df))
    gene_group[df$gene %in% c("Dnmt1","Eed","Ezh2","Kmt2a","Kmt2d")] = "chromatin modifiers"
    gene_group[df$gene %in% c("Ctcf","Gabpa","Nrf1")] = "ES TFs"
    gene_group[df$gene %in% c("NTC")] = "NTCs"
    df$gene_group = factor(gene_group, levels = c("NTCs","chromatin modifiers","ES TFs","developmental TFs"))
    
    df$comparing = factor(df$comparing, levels = c("ESC","EB"))
    
    p = ggplot(df, aes(gene_group, log2_fc, fill = gene_group)) + geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2) + facet_grid(.~comparing) +
        labs(x="", y="Log2 gRNA FC to plasmid library") +
        theme_classic(base_size = 10) +
        theme(legend.position="none") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values = pilot_gene_group_color_code) +
        theme(axis.text.x = element_text(color="black", angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(color="black")) 
    
    pdf(paste0(work_path, "/processing/pilot_screen_sgRNA/", run_id, "/LogFC_EB_vs_plasmid_boxplot.pdf"), 4, 5)
    print(p)
    dev.off()
}









