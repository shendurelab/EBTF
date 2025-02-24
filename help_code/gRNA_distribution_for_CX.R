#####Assessing gRNA distribution across cell types




#*********************************************************************************************************************************************
######Chi-square test for cluster-based gRNA distribution; note that several of these functions are from Jose


### Expectation maximiation model to correct for different efficiencies across sgRNAs

### Run on cds object with cluster information in metadata column called 'Clusters'


### load functions
get.guide.weights = function(mat, ntc.dist, n.iterations = 30) {
  n.guides = nrow(mat)
  n.cells = rowSums(mat)
  empirical.dist = sweep(mat, 1, n.cells, "/") #proportion in each cluster for each guide of a target
  
  lof.prop = rep(0.5, n.guides)
  expected.n.lof = n.cells * lof.prop
  
  for (i in 1:n.iterations) {
    lof.dist = sapply(1:n.guides, function(guide) {
      p = lof.prop[guide]
      (empirical.dist[guide,] - (1-p) * ntc.dist) / p    #(prop guide in cluster - 0.5*prop ntc in cluster)/0.5 -> = no adjustment if same as NTC; smaller fraction if less than (can be negative); larger if more than
    })
    
    lof.dist = rowSums(sweep(lof.dist, 2, expected.n.lof / sum(expected.n.lof), "*"))  #multiply with fraction of cells per guide (at 50% KO efficiency)
    lof.dist = ifelse(lof.dist < 0, 0, lof.dist)
    lof.dist = lof.dist / sum(lof.dist)  #sum over adjusted proportions for all guides
    
    lof.prop = sapply(1:n.guides, function(guide) {
      optimize(function(p) dmultinom(mat[guide,], prob = p * lof.dist + (1-p) * ntc.dist, log = T),
               c(0.0, 1.0), maximum = T)$maximum
    })
    
    expected.n.lof = n.cells * lof.prop   #different value for KO efficiency are tried; to explain enrichment/depletion values of sum of 3 guides
  }
  
  return(lof.prop)   #loss of function fraction for each of 3 guides
}


analysis.guides = 
  (colData(cds) %>% 
     as.data.frame() %>%
     filter(guide_count == 1) %>% 
     group_by(gene, barcode) %>%
     summarize(n.guide.cells = n()) %>% 
     group_by(gene) %>% 
     mutate(n.target.cells = sum(n.guide.cells)) %>%
     filter(n.guide.cells >= 10) %>% 
     ungroup())$barcode #3036; selecting guides with at least 10 cells


analysis.targets = as.data.frame(colData(cds) %>% 
                                   as.data.frame() %>%
                                   filter(gene != "NONTARGETING") %>%
                                   group_by(gene) %>% 
                                   summarize(n.cells = n(),
                                             n.guides = length(intersect(unique(barcode), analysis.guides))) %>%
                                   filter(n.cells >= 15, n.guides >= 1) %>% dplyr::select(gene))[,1] #1246; at least 15 cells per target
                                   
target.region.mat = reshape2::acast(
  colData(cds) %>%
    as.data.frame() %>%
    filter(barcode %in% analysis.guides | gene == "NONTARGETING") %>%
    mutate(dummy = 1) %>% dplyr::select(gene, Clusters, dummy),
  gene ~ Clusters, value.var = "dummy", fun.aggregate = sum, fill = 0) #simple count matrix of cells per target in each cluster

target.to.guide.map <- sapply(analysis.targets, function(x){
  
  sort(unique(as.data.frame(colData(cds) %>%
                              as.data.frame() %>%
                              filter(gene == x, barcode %in% analysis.guides) %>%
                              dplyr::select(barcode))[, 1]))
  
}) #list of targets with corresponding gRNAs (that meet cutoff criteria)


NTC.guides <- unique(colData(cds)[!is.na(colData(cds)$gene),][colData(cds)[!is.na(colData(cds)$gene),]$gene == "NONTARGETING",]$barcode)

guide.region.mat = reshape2::acast(
  colData(cds) %>% 
    as.data.frame() %>% 
    filter(barcode %in% analysis.guides[!(analysis.guides %in% NTC.guides)]) %>%
    #mutate(dummy = 1) %>% dplyr::select(barcode, dp_cluster, dummy),
    #barcode ~ dp_cluster, value.var = "dummy", fun.aggregate = sum, fill = 0)
    mutate(dummy = 1) %>% dplyr::select(barcode, Clusters, dummy),
  barcode ~ Clusters, value.var = "dummy", fun.aggregate = sum, fill = 0) #same count table as above just for guides, not targets

ntc.distribution = target.region.mat["NONTARGETING",]
ntc.distribution = ntc.distribution / sum(ntc.distribution) #fraction of NTCs across each cluster

weighted.target.region.mat = t(sapply(analysis.targets,
                                      function(target) {
                                        guides = target.to.guide.map[[target]]
                                        if (length(guides) == 1) {
                                          return(target.region.mat[target,])
                                        } else {
                                          mat = guide.region.mat[guides,]
                                          guide.weights = get.guide.weights(mat, ntc.distribution)
                                          guide.weights = guide.weights / max(guide.weights)
                                          
                                          print(target)
                                          print(round(guide.weights, 3))
                                          return(round(colSums(sweep(mat, 1, guide.weights, "*"))))
                                        }
                                      }))

NTC.region.mat <- matrix(target.region.mat[row.names(target.region.mat) == "NONTARGETING",],
                         nrow = 1)
row.names(NTC.region.mat) <- "NONTARGETING"

weighted.target.region.mat <- rbind(weighted.target.region.mat,
                                    NTC.region.mat)

NTC.region.p <- colData(cds) %>%
  as.data.frame() %>%
  group_by(gene, Clusters) %>%
  summarize(n = n()) %>%
  tidyr::complete(Clusters, fill = list(n = 0.1)) %>%
  filter(gene == "NONTARGETING")

ntc.distribution = weighted.target.region.mat["NONTARGETING",]    
ntc.distribution = ntc.distribution / sum(ntc.distribution)

set.seed(42)
    initial.target.level.chisq.pval = sapply(
      analysis.targets, function(target) {
        suppressWarnings({chisq.test(
          weighted.target.region.mat[target,],
          p = NTC.region.p$n,
          simulate.p.value = F, rescale.p = T, B = 1000)$p.value})
      })
    
initial.target.level.chisq.qval <- p.adjust(initial.target.level.chisq.pval, method = "BH") #this is the result of the ch-square test



#***********************************************************************************************************************************************





#Subsample NTCs so that we end up with similar amount of cells per barcode for targeting and nontargeting guides; otherwise NTCs don't serve as good control



