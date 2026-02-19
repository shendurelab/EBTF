
### The power to detect changes in cell type proportions upon perturbation of a target gene, relative to non-targeting controls, in a monoclonal organoid screen, is anticipated to be dependent on factors including:
### 
### 1) the abundance of the impacted cell type within each EB
### 2) the number of replicate ‘individuals’ (i.e., monoclonal EBs or clonotypes)
### 3) the effect size of the perturbation on the impacted cell type
### 4) the success of the perturbation (i.e. knockdown efficiency or penetrance)
### 5) the number of cells profiled per EB
### 6) the choice of statistical framework

### contact: cxqiu@uw.edu

########################################
### Method-1: Hooke based simulation ###
########################################

### Loading necessary function and gene meta data

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

library(tidyr)
library(forcats)
library(VGAM)
work_path = "Your_work_path"

### I used our monoclonal EBs (two conditions, pooled and arrayed) as examples, but this analysis can be extended to broader cases.
### The basic idea is to generate a table of cell counts across NTC samples and cell types using the real data.

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
### 102,120 cells, including 12 cell types

### "Arrayed" condition assignments (6493 cells, 19 NTC samples)
gRNA_assign_1 = readRDS(paste0(work_path, 
                               "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep1_used_for_hooke.rds"))
gRNA_assign_1 = gRNA_assign_1[[1]] %>% left_join(gRNA_assign_1[[2]], by = "clonotype_id") %>%
    filter(target == "NONTARGETING", !is.na(target)) %>% mutate(clonotype = paste0("arrayed_", clonotype_id))

### "Pooled" condition assignments (1337 cells, 5 NTC samples)
gRNA_assign_2 = readRDS(paste0(work_path, 
                               "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep2_used_for_hooke.rds"))
gRNA_assign_2 = gRNA_assign_2[[1]] %>% left_join(gRNA_assign_2[[2]], by = "clonotype_id", relationship = "many-to-many") %>%
    filter(target == "NONTARGETING", !is.na(target)) %>% mutate(clonotype = paste0("pooled_", clonotype_id))

gRNA_assign = rbind(gRNA_assign_1[,c("cell_id", "clonotype")], 
                    gRNA_assign_2[,c("cell_id", "clonotype")])

### two NTC are larger than mean + 2*sd, so excluding them
### 4690 cells, 22 NTC samples
clonotype_cell_num = gRNA_assign %>% group_by(clonotype) %>% tally()
clonotype_cell_num = clonotype_cell_num %>%
    filter(n < mean(clonotype_cell_num$n) + 2*sd(clonotype_cell_num$n))
gRNA_assign = gRNA_assign[gRNA_assign$clonotype %in% clonotype_cell_num$clonotype,]

coldata_df = pd %>% filter(cell_id %in% gRNA_assign$cell_id) %>%
    select(cell_id, celltype) %>% left_join(gRNA_assign, by = "cell_id")

types = coldata_df %>% group_by(celltype) %>% tally() %>% 
    arrange(desc(n)) %>% pull(celltype)

prob = coldata_df %>% group_by(celltype) %>% tally() %>%
    mutate(percent = (n/sum(n))) %>% arrange(-percent) %>%
    filter(celltype %in% types) %>% pull(percent)


##########################################################################################
### Applying the Hooke analysis to perform cell-type-composition test for each simulations

library(PLNmodels, lib.loc = "/net/gs/vol1/home/cxqiu/R/x86_64-pc-linux-gnu-library/4.3")
library(hooke)

obj = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned.rds"))
pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
count = GetAssayData(obj, slot = "counts")
keep = pd$cell_id %in% gRNA_assign$cell_id
pd_sub = pd[keep,]
count_sub = count[,keep]
pd_sub = pd_sub %>% select(cell_id, celltype) %>% left_join(gRNA_assign, by = "cell_id") %>% as.data.frame()
rownames(pd_sub) = as.vector(pd_sub$cell_id)
colnames(pd_sub) = c("cell_id", "cell_group", "embryo")

cds = new_cell_data_set(count_sub,
                        cell_metadata = pd_sub)

ccs = new_cell_count_set(cds, 
                         sample_group = "embryo", 
                         cell_group = "cell_group")

source("/help_code/power_analysis_hooke.R")
library(tibble)

effect_size_list = c(-0.75,-0.5,-0.25,-0.1,0.1,0.25,0.5,0.75)
num_embryos_list = c(10, 20, 30, 40, 50, 100, 200)
sample_size_list = c(50, 100, 200, 500, 1000)
success_ratio_list = c(0.5, 0.75, 1)

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1]) ### 1:10

for(sample_size in sample_size_list){
    result = NULL
    for(celltype_i in types){
        for(effect_i in effect_size_list){
            for(num_i in num_embryos_list){
                for(success_ratio_i in success_ratio_list){
                    print(paste0(sample_size, "/", celltype_i, "/", effect_i, "/", num_i))
                    res = run_comparison(ccs, 
                                         embryo_size = sample_size, 
                                         effect_size = (-1) * effect_i,
                                         num_embryos = num_i,
                                         cell_types = c(celltype_i),
                                         success_ratio = success_ratio_i,
                                         random.seed = kk, 
                                         num_threads = 4) 
                    result = rbind(result, data.frame(celltype = celltype_i,
                                                      effect_size = effect_i,
                                                      num_embryos = num_i,
                                                      success_ratio = success_ratio_i,
                                                      delta_log_abund = res$delta_log_abund[res$cell_group == celltype_i],
                                                      delta_p_value = res$delta_p_value[res$cell_group == celltype_i],
                                                      delta_q_value = res$delta_q_value[res$cell_group == celltype_i]))
                }
            }
        }
    }
    saveRDS(result, paste0(work_path, "/analysis/revision/power_analysis/iter_", kk, "/result_hooke_simulation_", sample_size, ".rds"))
}






####################################################
### Method-2: beta-binomial regression framework ###
####################################################

### Loading necessary function and gene meta data

source("help_script.R")
mouse_gene = read.table("mouse.GRCm38.p6.geneID.txt", header = T, as.is = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene_sub = mouse_gene[mouse_gene$chr %in% paste0("chr", c(1:19, "X", "Y", "M")),]

library(tidyr)
library(forcats)
library(VGAM)
work_path = "Your_work_path"

### I used our monoclonal EBs (two conditions, pooled and arrayed) as examples, but this analysis can be extended to broader cases.
### The basic idea is to generate a table of cell counts across NTC samples and cell types using the real data.

pd = readRDS(paste0(work_path, "/analysis/monoclonal_EB_TF_screen/obj_aligned_pd.rds"))
### 102,120 cells, including 12 cell types

### "Arrayed" condition assignments (6493 cells, 19 NTC samples)
gRNA_assign_1 = readRDS(paste0(work_path, 
                               "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep1_used_for_hooke.rds"))
gRNA_assign_1 = gRNA_assign_1[[1]] %>% left_join(gRNA_assign_1[[2]], by = "clonotype_id") %>%
    filter(target == "NONTARGETING", !is.na(target)) %>% mutate(clonotype = paste0("arrayed_", clonotype_id))

### "Pooled" condition assignments (1337 cells, 5 NTC samples)
gRNA_assign_2 = readRDS(paste0(work_path, 
                               "/analysis/monoclonal_EB_TF_screen/cell_to_mEB/cell_clonotype_gRNA_rep2_used_for_hooke.rds"))
gRNA_assign_2 = gRNA_assign_2[[1]] %>% left_join(gRNA_assign_2[[2]], by = "clonotype_id", relationship = "many-to-many") %>%
    filter(target == "NONTARGETING", !is.na(target)) %>% mutate(clonotype = paste0("pooled_", clonotype_id))

gRNA_assign = rbind(gRNA_assign_1[,c("cell_id", "clonotype")], 
                    gRNA_assign_2[,c("cell_id", "clonotype")])

### two NTC are larger than mean + 2*sd, so excluding them
### 4690 cells, 22 NTC samples
clonotype_cell_num = gRNA_assign %>% group_by(clonotype) %>% tally()
clonotype_cell_num = clonotype_cell_num %>%
    filter(n < mean(clonotype_cell_num$n) + 2*sd(clonotype_cell_num$n))
gRNA_assign = gRNA_assign[gRNA_assign$clonotype %in% clonotype_cell_num$clonotype,]

coldata_df = pd %>% filter(cell_id %in% gRNA_assign$cell_id) %>%
    select(cell_id, celltype) %>% left_join(gRNA_assign, by = "cell_id")

types = coldata_df %>% group_by(celltype) %>% tally() %>% 
    arrange(desc(n)) %>% pull(celltype)

prob = coldata_df %>% group_by(celltype) %>% tally() %>%
    mutate(percent = (n/sum(n))) %>% arrange(-percent) %>%
    filter(celltype %in% types) %>% pull(percent)


#########################################
### using beta binomial regression method

coldata_summary = coldata_df %>% 
    filter(celltype %in% types) %>% 
    group_by(clonotype, celltype) %>% tally() %>% rename(cells = n) %>%
    ungroup() %>%
    group_by(clonotype) %>%
    mutate(total_cells = sum(cells))

######################################################################
### step-2, fit a dirichlet model to simulate cell type abundances ###
######################################################################

### The Dirichlet model specifies probabilistically how many purchases each consumer 
### makes in a time-period and which brand is bought on each occasion. It combines 
### both purchase incidence and brand-choice aspects of buyer behaviour into one model.

coldata_fill = coldata_summary %>%
    tidyr::pivot_wider(id_cols = clonotype, names_from = celltype, values_from = cells, values_fill
                       = list(cells = 0)) %>%
    tidyr::gather(key = celltype, value = cells, -clonotype) %>%
    mutate(pseudo_cells = cells + 1) %>%
    group_by(clonotype) %>%
    mutate(total_cells = sum(pseudo_cells)) %>%
    mutate(cell_prop = pseudo_cells/total_cells) %>%
    ungroup()

# make a table of cell type proportions across individuals
cell_prop_wide = coldata_fill %>%
    select(clonotype, celltype, cell_prop) %>%
    pivot_wider(id_cols = clonotype, names_from = celltype, values_from = cell_prop)
prop_mat = as.matrix(cell_prop_wide[,-1])

### update the replicate size in each group to 5, 10, 15, ..., 50
prop_mat = as.matrix(cell_prop_wide[sample(nrow(cell_prop_wide), 5),-1])
prop_mat = prop_mat[,types]

# check that the rows all sum to 1
unique(rowSums(prop_mat)) == 1

# what do the porportions look like?
sort(colMeans(prop_mat))

# fit dirichlet
dfit <- VGAM::vglm(prop_mat ~ 1, data = as.data.frame(prop_mat), family = dirichlet, trace = FALSE)

######################################################################################
### step-3, simulate cell type abundances for our fake clonotype total count data. ###
######################################################################################

sims = c(1:7)
sims_list = c(2,4,6,8,10,20,40)
emb_num = nrow(prop_mat) 
#set.seed(111)
wt_dat = lapply(sims, function(x){
    dsim <- VGAM::simulate.vlm(dfit, nsim = sims_list[x])
    aaa <- array(unlist(dsim), c(emb_num*sims_list[x], ncol(fitted(dfit)), sims_list[x]))
    aaa = aaa[,,1]
})

wt_sim_props = melt(wt_dat)
colnames(wt_sim_props) = c("clonotype", "celltype", "cell_prop", "sim")
wt_sim_props$clonotype_unq = paste(wt_sim_props$sim, wt_sim_props$clonotype, sep='-')
print(wt_sim_props %>% group_by(celltype) %>% summarize(mean_cell_prop = mean(cell_prop)))

# how many cells per clonotype? 
cell_totals = coldata_summary %>% select(clonotype, total_cells) %>% unique()
size = round(mean(cell_totals$total_cells)) ### 213

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1]) ### 1:10

for(size in c(50, 100, 200, 500, 1000)){
    print(size)

    # get standard deviation from clonotype totals in our experiment
    my.sd = sd(cell_totals$total_cells)

    sim.total.cells = round(rnorm(1000, mean = size, sd = my.sd))
    if(size == 50){
        sim.total.cells = sim.total.cells[sim.total.cells >= 20]
    } else if (size == 100) {
        sim.total.cells = sim.total.cells[sim.total.cells >= 50]
    } else {
        sim.total.cells = sim.total.cells[sim.total.cells >= 100]
    }
    sim.total.cells = sample(sim.total.cells, length(unique(wt_sim_props$clonotype_unq)))
    wt_tot_df = data.frame(unique(wt_sim_props$clonotype_unq), sim.total.cells, stringsAsFactors = F)
    colnames(wt_tot_df) = c('clonotype_unq', 'total_cells')

    wt_sim_count_df = wt_sim_props %>%
        left_join(wt_tot_df, by='clonotype_unq') %>%
        mutate(cells = round(total_cells*cell_prop)) #%>%
    #group_by(clonotype_unq) %>%
    #mutate(total_cells = sum(cells)) #maybe adjust total cells to deal with rounding disc repencies

    #hist(wt_sim_count_df$total_cells) # check for normal distribution of total counts

    wt_sim_filt = wt_sim_count_df %>%
        left_join(wt_sim_count_df %>%
                      select(sim, clonotype_unq) %>%
                      distinct() %>%
                      group_by(sim) %>%
                      dplyr::mutate(sim_emb_count = dplyr::row_number()) %>%
                      ungroup() %>%
                      select(clonotype_unq, sim_emb_count), by = "clonotype_unq") %>%
        mutate(celltype = paste("celltype", celltype, sep = "_")) %>%
        select(sim_emb_count, celltype, sim, clonotype_unq, cells, total_cells) %>%
        group_by(sim) %>%
        add_count() %>%
        mutate(sim_emb_number = n/max(wt_sim_count_df$celltype), genotype = "wt") %>%
        select(-n) %>%
        ungroup() %>%
        arrange(clonotype_unq)

    print(sum(wt_sim_filt$cells < 0))
    ### considering set a different seed, to avoid negative value in the simulation

    ###################################################################################
    ### step-4, to look at our ability to detect different effect sizes, we'll first 
    ### need to split this data to use as controls vs. "mutants".
    ###################################################################################

    #set.seed(123)
    mut_dat = lapply(sims, function(x){
        dsim <- VGAM::simulate.vlm(dfit, nsim = sims_list[x])
        aaa <- array(unlist(dsim), c(emb_num*sims_list[x], ncol(fitted(dfit)), sims_list[x]))
        aaa = aaa[,,1]
    })

    mut_sim_props = melt(mut_dat)
    colnames(mut_sim_props) = c("clonotype", "celltype", "cell_prop", "sim")
    mut_sim_props$clonotype_unq = paste(mut_sim_props$sim, mut_sim_props$clonotype, sep='-')
    print(mut_sim_props %>% group_by(celltype) %>% summarize(mean_cell_prop = mean(cell_prop)))

    mut.total.cells = round(rnorm(1000, mean = size, sd = my.sd))
    if(size == 50){
        mut.total.cells = mut.total.cells[mut.total.cells >= 20]
    } else if (size == 100){
        mut.total.cells = mut.total.cells[mut.total.cells >= 50]
    } else {
        mut.total.cells = mut.total.cells[mut.total.cells >= 100]
    }
    mut.total.cells = sample(mut.total.cells, length(unique(mut_sim_props$clonotype_unq)))
    mut_tot_df = data.frame(unique(mut_sim_props$clonotype_unq), mut.total.cells, stringsAsFactors = F)
    colnames(mut_tot_df) = c('clonotype_unq', 'total_cells')

    mut_sim_count_df = mut_sim_props %>%
        left_join(mut_tot_df, by='clonotype_unq') %>%
        mutate(cells = round(total_cells*cell_prop))
    mut_sim_filt = mut_sim_count_df %>%
        left_join(mut_sim_count_df %>%
                      select(sim, clonotype_unq) %>%
                      distinct() %>%
                      group_by(sim) %>%
                      dplyr::mutate(sim_emb_count = dplyr::row_number()) %>%
                      ungroup() %>%
                      select(clonotype_unq, sim_emb_count), by = "clonotype_unq") %>%
        mutate(celltype = paste("celltype", celltype, sep = "_")) %>%
        select(sim_emb_count, celltype, sim, clonotype_unq, cells, total_cells) %>%
        group_by(sim) %>%
        add_count() %>%
        mutate(sim_emb_number = n/max(mut_sim_props$celltype),
               genotype = "mut") %>%
        select(-n) %>%
        ungroup() %>%
        arrange(clonotype_unq)

    print(sum(mut_sim_filt$cells < 0))
    ### considering set a different seed, to avoid negative value in the simulation

    saveRDS(list(wt_sim_filt, mut_sim_filt, types, prob),
            paste0(work_path, "/analysis/revision/power_analysis/iter_", kk, "/simulated_data_beta_binomial_", size, ".rds"))

    ###################################################
    ### step-5, performing beta binomial regression ###
    ###################################################

    binom_effect = function(wt_df, mut_df, which_type = 1, effect_vec = effects, success_ratio = success_ratio){
        test_res = lapply(effect_vec,
                          FUN = function(x) {
                              # filter for cell type
                              wt_df = wt_df %>% filter(celltype == which_type)
                              head(wt_df)
                              mut_df = mut_df %>% filter(celltype == which_type)
                              head(mut_df)
                              # adjust mutant data
                              #adj_df = mut_df %>% mutate(cells = round(cells * (1-x)))
                              adj_df = mut_df %>% mutate(cells = ifelse(runif(n()) < success_ratio, round(cells * (1 - x)), cells))
                              # combine and run test
                              comb_df = rbind(wt_df, adj_df) %>% 
                                  mutate(genotype = fct_relevel(genotype, c(unique(wt_df$genotype))))
                              
                              comb_df$cells = as.integer(comb_df$cells)
                              comb_df$total_cells = as.integer(comb_df$total_cells)
                              
                              count_df = cbind(comb_df$cells, comb_df$total_cells - comb_df$cells)
                              
                              fit = vglm(count_df ~ genotype, data = comb_df, family = betabinomial, trace = F)
                              
                              fit_df = tidy.vglm(fit)[3,]})
        test_res
    }


    binom_type = function(wt_df, mut_df, effect_vec = effects, type_vec = type_vec, success_ratio = success_ratio, which_sim = 1) {
        celltype_res = lapply(type_vec, function(x){
            wt_df = wt_df %>% dplyr::filter(sim == which_sim)
            mut_df = mut_df %>% dplyr::filter(sim == which_sim)
            res = binom_effect(wt_df = wt_df, mut_df = mut_df, which_type = x, effect_vec = effect_vec, success_ratio = success_ratio)
            names(res) = effects
            res.bind = do.call(rbind, res)
            res.bind$effect_size = row.names(res.bind)
            res.bind$celltype = x
            res.bind
        })
        celltype_res
    }


    tidy.vglm = function(x, conf.int=FALSE, conf.level=0.95) {
        co <- as.data.frame(coef(summary(x)))
        names(co) <- c("estimate","std.error","statistic","p.value") 
        if (conf.int) {
            qq <- qnorm((1+conf.level)/2)
            co <- transform(co,
                            conf.low=estimate-qq*std.error,
                            conf.high=estimate+qq*std.error)
        }
        co <- data.frame(term=rownames(co),co) 
        rownames(co) <- NULL
        return(co)
    }


    compare_abundance = function(cell_df, wt_df){
        comb_df = rbind(as.data.frame(cell_df), as.data.frame(wt_df)) %>%
            mutate(genotype = fct_relevel(genotype, c(unique(wt_df$genotype))))
        cell.types = unique(comb_df$celltype)
        test_res = sapply(cell.types,
                          FUN = function(x) {
                              type_df = comb_df %>% filter(celltype == x)
                              count_df = cbind(type_df$cells, type_df$total_cells - type_df$cells)
                              fit =  VGAM::vglm(count_df ~ genotype, data = type_df, family = "betabinomial", trace = F)
                              fit_df = tidy.vglm(fit)[3,]}, USE.NAMES = T, simplify = F)
        test_res = do.call(rbind, test_res)
        test_res = test_res %>% tibble::rownames_to_column(var = "cell_group")
        test_res %>% arrange(desc(estimate))
    }

    ################################################################################
    ### step-5, Test multiple effect sizes with increasing numbers of clonotypes ###
    ################################################################################
    ### Using the functions defined above, we’ll now look at whether we can detect 
    ### signficant differences in cell abundance across multiple effect sizes, with 
    ### different numbers of clonotypes as input.

    emb_interval = 5
    embs = sims_list*emb_interval

    v1 = unlist(lapply(embs, function(x) seq(1, x, 1)))
    v2 = unlist(lapply(sims, function(x) rep(x, sims_list[x]*emb_interval)))

    emb_set = tibble(v2, v1) %>%
        mutate(emb_name = paste(v2, v1, sep = "-")) %>%
        pull(emb_name)

    # filter simulated datasets
    wt_sim = wt_sim_filt %>%
        filter(clonotype_unq %in% emb_set)

    mut_sim = mut_sim_filt %>%
        filter(clonotype_unq %in% emb_set)

    # effect sizes
    effects = c(0.1,0.25,0.5,0.75)
    success_ratio_list = c(0.5, 0.75, 1)

    # simulations and cell types
    final.res.bind.comb = NULL
    for(success_ratio in success_ratio_list){
        type_vec = unique(wt_sim$celltype)
        # calculate significant differences
        final.res = pbapply::pblapply(sims, function(x){
        res = binom_type(wt_df = wt_sim, mut_df = mut_sim,
                         type_vec = type_vec, success_ratio = success_ratio, which_sim = x)
        names(res) = sims
        res.bind = do.call(rbind, res)
        res.bind$sim = x
        res.bind
        })
        final.res.bind = do.call(rbind, final.res)
        head(final.res.bind)
        final.res.bind$emb_num = sims_list[final.res.bind$sim] * emb_interval
        final.res.bind$success_ratio = success_ratio
        final.res.bind.comb = rbind(final.res.bind.comb, final.res.bind)
    }

    saveRDS(final.res.bind.comb, paste0(work_path, "/analysis/revision/power_analysis/iter_", kk, "/result_beta_binomial_simulation_", size, ".rds"))


}






