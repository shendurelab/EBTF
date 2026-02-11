
#' returns a cell group x embryo count matrix 
#' @param ccs
#' @return a cell group by embryo count matrix
get_cell_count_wide <- function(ccs) {
  # cell group x embryo
  cell_count_wide = counts(ccs) %>%
    as.matrix %>%
    as.data.frame
  
  return(cell_count_wide)
}

#' takes in a ccs and outputs a matrix of
#' cell_type x embryo of proportions
#' @param cell_count_wide a cell group x embryo count matrix
#' @param pseudocount
get_prop_mat <- function(cell_count_wide,
                         pseudocount = 1) {
  
  prop_mat = cell_count_wide %>%
    rownames_to_column("cell_group") %>%
    tidyr::pivot_longer(-cell_group, names_to = "embryo", values_to = "count") %>%
    mutate(count = count + pseudocount) %>%
    group_by(embryo) %>%
    mutate(prop = (count)/sum(count)) %>%
    select(-count) %>%
    tidyr::pivot_wider(names_from = cell_group, values_from = prop) %>%
    tibble::column_to_rownames("embryo") %>%
    as.matrix
  
  return(prop_mat)
  
}

#' fits a drichlet model based off proportions
#' @param prop_mat
#' @param model_formula_str
#' @param trace
fit_drichlet <- function(prop_mat,
                         model_formula_str = "~ 1",
                         trace = FALSE) {
  
  model_formula_str = paste("prop_mat", model_formula_str)
  model_formula = as.formula(model_formula_str)
  
  dfit <- VGAM::vglm(model_formula,
                     data = as.data.frame(prop_mat),
                     family = "dirichlet",
                     trace = trace)
  
  return(dfit)
}


#' simulates a fake ccs based off wt data
#' @param ccs 
#' @param embryo_size how many cells per embryo
#' @param num_embryos how many embryo replicates
#' @param genotype
#' @param random.seed 
simulate_ccs = function(ccs, 
                        embryo_size = 1000,
                        random.seed = 111, 
                        genotype = "WT", 
                        num_embryos = NULL) {
  
  cell_count_wide = get_cell_count_wide(ccs)
  prop_mat = get_prop_mat(cell_count_wide)
  prop_mat = prop_mat[sample(nrow(prop_mat), 5),]
  dfit = fit_drichlet(prop_mat)
  
  if (is.null(num_embryos)) {
    num_embryos = ncol(cell_count_wide)
  }
  num_cell_groups = nrow(cell_count_wide)
  
  num_sims = ceiling(num_embryos / nrow(prop_mat))
  
  dsim <- VGAM::simulate.vlm(dfit, nsim = num_sims, random.seed = random.seed)
  # this returns a (num_embryos x num_cell_groups) x num_sims matrix
  # make into a num_embryos x num_cell_groups x num_sims matrix 
  aaa <- array(unlist(dsim), c(num_embryos, num_cell_groups, num_sims))
  
  aaa_list = lapply(1:num_sims, function(i){
    # slice by dims
    aaa[,,i]
  })
  
  bbb = do.call(rbind, aaa_list) 
  bbb = bbb[1:num_embryos,]
  print(apply(bbb, 2, mean))
  
  # how many total cells does each embryo have
  cell_totals = colSums(cell_count_wide)
  # get standard deviation from embryo totals in our experiment
  my.sd = sd(cell_totals)
  
  # simulate total cells
  # sim.total.cells = round(rnorm(num_embryos, mean = embryo_size, sd = my.sd))
  # sim.total.cells = rpois(n = num_embryos, lambda = embryo_size)
  sim.total.cells = round(rnorm(1000, mean = embryo_size, sd = my.sd))
  sim.total.cells = sim.total.cells[sim.total.cells >= 100]
  sim.total.cells = sample(sim.total.cells, num_embryos)
  
  sim_props = reshape2::melt(bbb)
  colnames(sim_props) = c("embryo", "cell_group_number", "cell_prop")
  cell_group_names = rownames(cell_count_wide)
  df = data.frame("cell_group_number" = 1:length(cell_group_names),
                  "cell_group" = cell_group_names)
  sim_props = sim_props %>% left_join(df, by = "cell_group_number")
  
  tot_df = data.frame(unique(sim_props$embryo), sim.total.cells, stringsAsFactors = F)
  colnames(tot_df) = c('embryo', 'total_cells')
  
  
  # convert to counts
  sim_count_df = sim_props %>%
    # left_join(tot_df, by='embryo_unq') %>%
    left_join(tot_df, by='embryo') %>%
    mutate(count = round(total_cells*cell_prop)) %>%
    select(-c(total_cells, cell_prop, cell_group_number))
  
  sim_count_wide = sim_count_df %>%
    pivot_wider(names_from = "embryo", values_from = count) %>%
    tibble::column_to_rownames("cell_group") %>%
    as.matrix
  
  colnames(sim_count_wide) = paste0("embryo-", genotype, "-",colnames(sim_count_wide))
  # make a cell_count_cds
  covariates_df = data.frame("embryo" = colnames(sim_count_wide), 
                             genotype = genotype, 
                             effect = 0, 
                             embryo_size = embryo_size,
                             random.seed = random.seed, 
                             mutated = F)
  rownames(covariates_df) = covariates_df$embryo
  cell_count_cds = monocle3::new_cell_data_set(sim_count_wide,
                                               cell_metadata=covariates_df)
  sim_ccs = suppressMessages(
    methods::new("cell_count_set",
                 cell_count_cds,
                 cds = new_cell_data_set(hooke:::empty_sparse_matrix(format="C")), 
                 cds_coldata = tibble(), 
                 cds_reduced_dims = SimpleList(),
                 info = SimpleList())
  )
  
  return(sim_ccs)
  
}

#' mutates a ccs object
#'@param sim_ccs a simulated ccs object
#'@param cell_types a list of cell types to mutate, if null will randomly sample
#'@param num_cell_types integer how many cell types to randomly sample
#'@param effect_size how much to mutate cell types by
mutate_ccs = function(sim_ccs, 
                      cell_types = NULL,
                      num_cell_types = 10,
                      effect_size = 1.0, 
                      success_ratio = 0.8,
                      random.seed = 111) {
  
  sim_count_wide = get_cell_count_wide(sim_ccs)
  
  if (is.null(cell_types)) {
    set.seed(random.seed)
    all_cell_types = rownames(sim_count_wide)
    cell_types = sample(all_cell_types, num_cell_types)
  }
  
  sim_count_df = sim_count_wide %>% 
                 as.data.frame %>% 
                 rownames_to_column("cell_group") %>%
                 pivot_longer(-cell_group,
                             names_to = "embryo",
                             values_to = "count")

  sim_count_df = sim_count_df %>%
  mutate(count = ifelse(
    cell_group %in% cell_types & runif(n()) < success_ratio,
    round(count * (1 - effect_size)),
    count
  ))
  
  sim_count_wide = sim_count_df %>%
    pivot_wider(names_from = "embryo",
                values_from = count) %>%
    tibble::column_to_rownames("cell_group") %>%
    as.matrix()
  
  covariates_df = colData(sim_ccs)
  
  covariates_df$effect_size = effect_size
  covariates_df$mutated_cell_types = list(cell_types)
  rownames(covariates_df) = colnames(sim_ccs) 
  colnames(sim_count_wide) = colnames(sim_ccs) 
  
  cell_count_cds = monocle3::new_cell_data_set(sim_count_wide,
                                               cell_metadata=covariates_df)
  
  # now turn it back into a cell count ccs
  mt_ccs = suppressMessages(
            methods::new("cell_count_set",
                        cell_count_cds,
                        cds = new_cell_data_set(hooke:::empty_sparse_matrix(format="C")), 
                        cds_coldata = tibble(), 
                        cds_reduced_dims = SimpleList(),
                        info = SimpleList())
  )
  return(mt_ccs)
  
}

#' combine two ccs objects
#' @param ccs_list 
combine_ccs = function(ccs_list) {
  
  comb_counts_list = lapply(ccs_list, function(ccs) {
    counts(ccs)
  })
  comb_counts = do.call(cbind, comb_counts_list)
  
  comb_coldata_list = lapply(ccs_list, function(ccs) {
    colData(ccs) %>% as.data.frame
  })
  
  comb_coldata = bind_rows(comb_coldata_list)
  
  cell_count_cds = monocle3::new_cell_data_set(comb_counts,
                                               cell_metadata=comb_coldata)
  
  comb_ccs = suppressMessages(
    methods::new("cell_count_set",
                 cell_count_cds,
                 cds = new_cell_data_set(hooke:::empty_sparse_matrix(format="C")), 
                 cds_coldata = tibble(), 
                 cds_reduced_dims = SimpleList(),
                 info = SimpleList())
  )
  return(comb_ccs)
}

#' @param ccs cell_count_set made from real data
#' @param embryo_size
#' @param random.seed
#' @param cell_types list of cell types to apply mutation to 
#' @param effect_size defaults to 1 
run_comparison = function(ccs, 
                          embryo_size = 1000, 
                          effect_size = 1,
                          num_embryos = 8,
                          cell_types = c(),
                          success_ratio = 0.8,
                          random.seed = 111, 
                          num_threads = 4) {
  
  wt_ccs = simulate_ccs(ccs,
                        embryo_size = embryo_size,
                        num_embryos = num_embryos,
                        random.seed = random.seed, genotype = "WT")

  mt_ccs = simulate_ccs(ccs,
                        embryo_size = embryo_size,
                        num_embryos = num_embryos,
                        random.seed = random.seed*2, genotype = "MT")

  mt_ccs = mutate_ccs(mt_ccs, cell_types, effect_size = effect_size, success_ratio = success_ratio)

  ccs = combine_ccs(ccs_list = list(wt_ccs, mt_ccs))
  ccm = new_cell_count_model(ccs, 
                             main_model_formula_str = "~ genotype", 
                             penalize_by_distance = F, 
                             num_threads = num_threads, 
                             vhat_method = "bootstrap")
  
  cond_wt = estimate_abundances(ccm, tibble(genotype = "WT"))
  cond_mt = estimate_abundances(ccm, tibble(genotype = "MT"))
  cond_wt_v_mt_tbl = compare_abundances(ccm, cond_wt, cond_mt)
  return(cond_wt_v_mt_tbl)
}

