#' Observed beta diversity 
#' 
#' \code{obsbeta3D}: Observed beta diversity with order q
#' 
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a \code{matrix/data.frame} (species by assemblages), or a \code{list} of \code{matrices/data.frames}, each matrix represents species-by-assemblages abundance matrix; see Note 1 for examples.\cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a \code{list} of \code{matrices/data.frames}, each matrix represents species-by-sampling units; see Note 2 for an example.
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is c(0, 1, 2).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being \code{0} (non-detection) or \code{1} (detection).
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter \code{0} to skip the bootstrap procedures. Default is \code{20}. If more accurate results are required, set \code{nboot = 100} (or \code{200}).
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{0.95}.
#' @param PDtree (required only when \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage. 
#' @param PDreftime (required only when \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of PDtree).  
#' @param PDtype (required only when \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{"meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param FDdistM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled assemblage. 
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_value"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}.  
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_value"}), a numerical vector between 0 and 1 specifying tau value (threshold level). If \code{NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled assemblage (i.e., quadratic entropy). 
#' @param FDcut_number (required only when \code{diversity = "FD"} and \code{FDtype = "AUC"}), a numeric number to split zero to one into several equal-spaced length. Default is 30.
#' 
#' @import tidyverse
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import ape
#' @import phytools
#' @import phyclust
#' @import tidytree
#' @import colorRamps
#' @import iNEXT.3D
#' @import future.apply
#' @import ade4
#' @import tidyr
#' @import tibble
#' 
#' @return A list of seven lists with three-diversity and four-dissimilarity.
#' 
#' @examples
#' ## Taxonomic diversity for abundance data
#' data(beetle_abu)
#' output1 = obsbeta3D(data = beetle_abu, diversity = 'TD', datatype = 'abundance', nboot = 20, conf = 0.95)
#' output1
#' 
#' 
#' ## Taxonomic diversity for incidence data
#' data(beetle_inc)
#' output2 = obsbeta3D(data = beetle_inc, diversity = 'TD', datatype = 'incidence_raw', nboot = 20, conf = 0.95)
#' output2
#' 
#' 
#' ## Phylogenetic diversity for abundance data
#' data(beetle_abu)
#' data(beetle_tree)
#' output3 = obsbeta3D(data = beetle_abu, diversity = 'PD', datatype = 'abundance', nboot = 20, conf = 0.95, PDtree = beetle_tree, PDreftime = NULL, PDtype = 'PD')
#' output3
#' 
#' 
#' ## Phylogenetic diversity for incidence data
#' data(beetle_inc)
#' data(beetle_tree)
#' output4 = obsbeta3D(data = beetle_inc, diversity = 'PD', datatype = 'incidence_raw', nboot = 20, conf = 0.95, PDtree = beetle_tree, PDreftime = NULL, PDtype = 'PD')
#' output4
#' 
#' 
#' ## Functional diversity for abundance data under single threshold
#' data(beetle_abu)
#' data(beetle_distM)
#' output5 = obsbeta3D(data = beetle_abu, diversity = 'FD', datatype = 'abundance', nboot = 20, conf = 0.95, FDdistM = beetle_distM, FDtype = 'tau_value', FDtau = NULL)
#' output5
#' 
#' 
#' ## Functional diversity for incidence data under single threshold
#' data(beetle_inc)
#' data(beetle_distM)
#' output6 = obsbeta3D(data = beetle_inc, diversity = 'FD', datatype = 'incidence_raw', nboot = 20, conf = 0.95, FDdistM = beetle_distM, FDtype = 'tau_value', FDtau = NULL)
#' output6
#' 
#' 
#' ## Functional diversity for abundance data with thresholds integrating from 0 to 1
#' data(beetle_abu)
#' data(beetle_distM)
#' output7 = obsbeta3D(data = beetle_abu, diversity = 'FD', datatype = 'abundance', nboot = 10, conf = 0.95, FDdistM = beetle_distM, FDtype = 'AUC', FDcut_number = 30)
#' output7
#' 
#' 
#' ## Functional diversity for incidence data with thresholds integrating from 0 to 1
#' data(beetle_inc)
#' data(beetle_distM)
#' output8 = obsbeta3D(data = beetle_inc, diversity = 'FD', datatype = 'incidence_raw', nboot = 10, conf = 0.95, FDdistM = beetle_distM, FDtype = 'AUC', FDcut_number = 30)
#' output8
#' 
#' @export
obsbeta3D = function(data, diversity = 'TD', q = seq(0, 2, 0.25), datatype = 'abundance',
                     nboot = 20, conf = 0.95, PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD',
                     FDdistM = NULL, FDtype = 'AUC', FDtau = NULL, FDcut_number = 50) {
  if (datatype == 'abundance') {
    
    if( inherits(data, "data.frame") | inherits(data, "matrix") ) data = list(Region_1 = data)
    
    if(class(data) == "list"){
      if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
      Ns = sapply(data, ncol)
      data_list = data
    }
    
  }
  
  if (datatype == 'incidence_raw') {
    
    if(is.null(names(data))) region_names = paste0("Region_", 1:length(data)) else region_names = names(data)
    Ns = sapply(data, length)
    data_list = data
    
  }
  
  
  if (is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)
  
  if (diversity == 'FD' & FDtype == 'tau_value' & is.null(FDtau) == T) {
    if (datatype == 'abundance') {
      pdata <- sapply(data_list, rowSums) %>% rowSums
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
    } else if (datatype == 'incidence_raw') {
      pdata <- sapply(data_list, function(x) {tmp = Reduce('+', x); tmp[tmp > 1] = 1; rowSums(tmp) }) %>% rowSums
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
    }
    FDtau <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
  }
  
  if (diversity == 'PD') {
    
    if (datatype == 'abundance') pool.data = do.call(cbind, data_list) %>% rowSums
    
    if (datatype == 'incidence_raw') pool.data = do.call(cbind,lapply(data_list, function(x) do.call(cbind,x)) ) %>% rowSums
    
    pool.name = names(pool.data[pool.data>0])
    tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
    mytree = drop.tip(PDtree, tip)
    H_max = get.rooted.tree.height(mytree)
    
    if(is.null(PDreftime)) { reft = H_max
    } else if (PDreftime <= 0) { stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
    } else { reft = PDreftime }
    
  }
  
  for_each_region = function(data, region_name, N) {
    
    #data
    if (datatype == 'abundance') {
      
      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector
      
    }
    
    if (datatype == 'incidence_raw') {
      
      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all regions must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)
      
      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))
      
      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)
      
      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
      
    }
    
    
    
    if (diversity == 'TD') {
      
      if (datatype == 'abundance') {
        
        gamma = AO3D(as.numeric(data_gamma), 'TD', q = q, datatype = "abundance", nboot = 0, method = 'Observed')
        
        alpha = AO3D(as.numeric(data_alpha), 'TD', q = q, datatype = "abundance", nboot = 0, method = 'Observed')
        
      }
      
      if (datatype == 'incidence_raw') {
        
        gamma = AO3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0, method = 'Observed')
        
        alpha = AO3D(as.numeric(data_alpha_freq), diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0, method = 'Observed')
        
        
      }
      
      gamma = gamma[,c(2,1)] %>% set_colnames(c('Diversity', 'Order.q'))
      
      alpha = alpha[,c(2,1)] %>% set_colnames(c('Diversity', 'Order.q'))
      
      alpha$Diversity = alpha$Diversity / N
      
      beta = alpha
      beta$Diversity = gamma$Diversity/alpha$Diversity
      C = beta %>% mutate(Diversity = ifelse(Order.q==1, log(Diversity)/log(N), (Diversity^(1-Order.q) - 1)/(N^(1-Order.q)-1)))
      U = beta %>% mutate(Diversity = ifelse(Order.q==1, log(Diversity)/log(N), (Diversity^(Order.q-1) - 1)/(N^(Order.q-1)-1)))
      V = beta %>% mutate(Diversity = (Diversity-1)/(N-1))
      S = beta %>% mutate(Diversity = (1/Diversity-1)/(1/N-1))
      
      if(nboot>1){
        
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))
            
            bootstrap_data_gamma = rowSums(bootstrap_sample)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 0]
            bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 0]
            
            gamma = AO3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", nboot = 0, method = 'Observed')
            
            alpha = AO3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", nboot = 0, method = 'Observed')
            
          }
          
          if (datatype == 'incidence_raw') {
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            
            raw = lapply(1:ncol(bootstrap_population), function(j){
              
              lapply(1:nrow(bootstrap_population), function(i) rbinom(n = n, size = 1, prob = bootstrap_population[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            
            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
            
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
            
            gamma = AO3D(bootstrap_data_gamma_freq, diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0, method = 'Observed')
            
            alpha = AO3D(bootstrap_data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", nboot = 0, method = 'Observed')
            
          }
          
          gamma = gamma$qD
          
          alpha = alpha$qD
          alpha = alpha / N
          
          beta = data.frame(Diversity = gamma/alpha, q)
          
          C = (beta %>% mutate(Diversity = ifelse(q == 1,log(Diversity)/log(N),(Diversity^(1 - q) - 1)/(N^(1 - q) - 1))))$Diversity
          U = (beta %>% mutate(Diversity = ifelse(q == 1,log(Diversity)/log(N),(Diversity^(q - 1) - 1)/(N^(q - 1) - 1))))$Diversity
          V = (beta %>% mutate(Diversity = (Diversity - 1)/(N - 1)))$Diversity
          S = (beta %>% mutate(Diversity = (1/Diversity - 1)/(1/N - 1)))$Diversity
          
          beta = beta$Diversity
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
        }) %>% abind(along = 3) %>% apply(1:2, sd)
        
      } else {
        
        se = matrix(0, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    if (diversity == 'PD') {
      
      if (datatype == 'abundance') {
        
        aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        
        gamma = iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, reft = reft, cal = "PD", nt = n) %>% 
          matrix(., ncol = 1) %>% data.frame() %>% cbind(q) %>% set_colnames(c('Diversity', 'Order.q'))
        
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[data[,i]>0,i]
          names(x) = rownames(data)[data[,i]>0]
          
          aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        
        qPDm = iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, reft = reft, cal = "PD", nt = n) %>%
          matrix(., ncol = 1)
        qPDm = qPDm/N
        alpha = qPDm %>% data.frame() %>% cbind(q) %>% set_colnames(c('Diversity', 'Order.q'))
        
      }
      
      if (datatype == 'incidence_raw') {
        
        aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        gamma = iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis=as.matrix(aL_table_gamma$branch.length), q = q, reft = reft, cal = "PD", nt = n) %>% 
          matrix(., ncol = 1) %>% data.frame() %>% cbind(q) %>% set_colnames(c('Diversity', 'Order.q'))
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[[i]]
          
          aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        alpha = (iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, reft = reft, cal = "PD", nt = n)/N) %>% 
          matrix(., ncol = 1) %>% cbind(q) %>% data.frame() %>% set_colnames(c('Diversity', 'Order.q'))
        
        
      }
      
      if (PDtype == 'meanPD') {
        gamma$Diversity = gamma$Diversity/reft
        alpha$Diversity = alpha$Diversity/reft
      }
      
      beta = alpha
      beta$Diversity = gamma$Diversity/alpha$Diversity
      
      C = beta %>% mutate(Diversity = ifelse(Order.q == 1, log(Diversity)/log(N), (Diversity^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      U = beta %>% mutate(Diversity = ifelse(Order.q == 1, log(Diversity)/log(N), (Diversity^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      V = beta %>% mutate(Diversity = (Diversity - 1)/(N - 1))
      S = beta %>% mutate(Diversity = (1/Diversity - 1)/(1/N - 1))
      
      if(nboot>1){
        
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            tree_bt = PDtree
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol = ncol(data))
            
            if ( nrow(p_bt) > nrow(data) & sum(unseen_p) > 0 ){
              
              unseen = unseen_p[which(rowSums(unseen_p) > 0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data), unseen_name)
              
              bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              x_bt = bootstrap_sample
              
              rownames(x_bt) = rownames(p_bt)
              
              if ( sum(x_bt[-(1:nrow(data)),])>0 ){
                
                g0_hat = apply(data, 2, function(x){
                  
                  n = sum(x)
                  f1 = sum(x == 1)
                  f2 = sum(x == 2)
                  
                  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
                  
                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  g1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  g0_hat = ifelse( g2 > ((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
                  if(is.na(g0_hat)) {g0_hat <- 0 }
                  g0_hat
                  
                })
                
                te = (x_bt[1:nrow(data),]*(data == 0))>0
                used_length = sapply(1:ncol(data), function(i) { 
                  
                  if (sum(te[,i]) == 0) return(0) else {
                    
                    iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
                    
                  }
                  
                })
                
                g0_hat = g0_hat - used_length
                g0_hat[g0_hat < 0] = 0
                
                unseen_sample = x_bt[-(1:nrow(data)),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = ncol(x_bt), byrow = T)
                
                L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (g0_hat[i] / nrow(unseen)) else 0 )
                
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0
                
                for (i in 1:length(L0_hat)){
                  
                  tip = list(edge = matrix(c(2,1),1,2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"
                  
                  tree_bt = tree_bt + tip
                  
                }
                
              } else {
                
                x_bt = x_bt[1:nrow(data),]
                p_bt = p_bt[1:nrow(data),]
                
              }
              
            } else {
              
              p_bt = p_bt[1:nrow(data),]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              rownames(x_bt) = rownames(data)
              
            }
            
            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]
            
            aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            gamma = as.vector(iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, reft = reft, cal = "PD", nt = n))
            
            
            aL_table_alpha = c()
            
            for (i in 1:N){
              
              x = x_bt[,i]
              names(x) = rownames(p_bt)
              x = x[x_bt[,i]>0]
              
              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            alpha = as.vector((iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, reft = reft, cal = "PD", nt = n)/N))
            
          }
          
          if (datatype == 'incidence_raw') {
            
            tree_bt = PDtree
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])),] %>% matrix(ncol=N)
            
            if ( nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p)>0 ){
              
              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data[[1]])),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data[[1]]), unseen_name)
              
              raw = lapply(1:ncol(p_bt), function(j){
                
                lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)
                
              })
              
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              
              if ( lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])),]) %>% do.call(sum,.)>0 ){
                
                R0_hat = sapply(data, function(x){
                  
                  nT = ncol(x)
                  Q1 = sum(rowSums(x)==1)
                  Q2 = sum(rowSums(x)==2)
                  
                  aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
                  
                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  R1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  R0_hat = ifelse( R2>((R1*Q2)/(2*Q1)) , ((nT-1)/nT)*(R1^2/(2*R2)) , ((nT-1)/nT)*(R1*(Q1-1)/(2*(Q2+1))) )
                  if(is.na(R0_hat)) { R0_hat <- 0 }
                  R0_hat
                  
                })
                
                te = (sapply(raw, rowSums)[1:nrow(data[[1]]),]*(sapply(data, rowSums) == 0)) > 0
                used_length = sapply(1:N, function(i) {
                  
                  if (sum(te[,i]) == 0) return(0) else {
                    
                    iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = raw[[i]][1:nrow(data[[1]]),], datatype = "incidence_raw", rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
                    
                  }
                  
                })
                
                R0_hat = R0_hat - used_length
                R0_hat[R0_hat < 0] = 0
                
                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = N, byrow = T)
                
                L0_hat = sapply(1:length(R0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (R0_hat[i] / nrow(unseen)) else 0 )
                
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0
                
                for (i in 1:length(L0_hat)){
                  
                  tip = list(edge = matrix(c(2,1), 1, 2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"
                  
                  tree_bt = tree_bt + tip
                  
                }
                
              } else raw = lapply(raw, function(i) i[1:nrow(data[[1]]),])
              
            } else {
              
              p_bt = p_bt[1:nrow(data[[1]]),]
              raw = lapply(1:ncol(p_bt), function(j){
                
                lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
                
              })
              
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              
            }
            
            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            
            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
            
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
            
            aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = bootstrap_data_gamma_raw, datatype = "incidence_raw", rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            gamma = as.vector(iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, reft = reft, cal = "PD", nt = n))
            
            
            aL_table_alpha = c()
            
            for (i in 1:N){
              
              x = raw[[i]]
              
              aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            alpha = as.vector((iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, reft = reft, cal = "PD", nt = n)/N))
            
          }
          
          if (PDtype == 'meanPD') {
            gamma = gamma/reft
            alpha = alpha/reft
          } 
          
          beta = data.frame(Diversity = gamma/alpha, Order.q = q)
          
          C = (beta %>% mutate(Diversity = ifelse(Order.q == 1, log(Diversity)/log(N), (Diversity^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Diversity
          U = (beta %>% mutate(Diversity = ifelse(Order.q == 1, log(Diversity)/log(N), (Diversity^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Diversity
          V = (beta %>% mutate(Diversity = (Diversity - 1)/(N - 1)))$Diversity
          S = (beta %>% mutate(Diversity = (1/Diversity - 1)/(1/N - 1)))$Diversity
          
          beta = beta$Diversity
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
        }) %>% abind(along = 3) %>% apply(1:2, sd)
        
      } else {
        
        se = matrix(0, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    if (diversity == 'FD') {
      
      FDdistM = as.matrix(FDdistM)
      
      FD_by_tau = function(data, distM, tau, datatype) {
        
        if (datatype == 'abundance') {
          
          zik = data
          zik = zik[rowSums(data)>0,]
          
          dij = distM
          dij = dij[rowSums(data)>0, rowSums(data)>0]
          
          dij[which(dij>tau, arr.ind = T)] = tau
          aik = (1 - dij/tau) %*% as.matrix(zik)
          positive_id = rowSums(aik)>0
          
          
          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_v = gamma_x/gamma_a
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          gamma = iNEXT.3D:::FD_mle(ai_vi_gamma, q) %>% as.vector
          
          
          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          
          alpha_v = alpha_x/alpha_a
          alpha_v = rep(gamma_v,N)
          
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          alpha = (iNEXT.3D:::FD_mle(ai_vi_alpha, q)/N) %>% as.vector
          
        }
        
        if (datatype == 'incidence_raw') {
          
          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D
          
          gamma_Y = data_gamma_freq[-1]
          
          dij = distM
          dij = dij[gamma_Y > 0, gamma_Y > 0]
          gamma_Y = gamma_Y[gamma_Y > 0]
          
          dij[which(dij>tau, arr.ind = T)] = tau
          gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          gamma_a[gamma_a > n] = n
          
          gamma_v = gamma_Y/gamma_a
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          gamma = iNEXT.3D:::FD_mle(ai_vi_gamma, q) %>% as.vector
          
          alpha_Y = data_2D[-1,]
          
          dij = distM
          dij = dij[rowSums(data_2D[-1,]) > 0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]
          
          dij[which(dij>tau, arr.ind = T)] = tau
          alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)
          
          alpha_a[alpha_a > n] = n
          alpha_a = as.vector(alpha_a)
          
          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a > 0]
          alpha_a = alpha_a[alpha_a > 0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          alpha = (iNEXT.3D:::FD_mle(ai_vi_alpha, q)/N) %>% as.vector
          
        }
        
        return(data.frame(gamma,alpha))
        
      }
      
      if (FDtype == 'tau_value'){
        
        if (datatype == 'abundance') {
          
          output = FD_by_tau(data, FDdistM, FDtau, datatype = 'abundance')
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(gamma, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
          alpha = data.frame(alpha, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
          beta = alpha
          beta$Diversity = gamma$Diversity/alpha$Diversity
          
        }
        
        if (datatype == 'incidence_raw') {
          
          output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, datatype='incidence_raw')
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(gamma, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
          alpha = data.frame(alpha, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
          beta = alpha
          beta$Diversity = gamma$Diversity/alpha$Diversity
          
        }
        
      }
      
      if (FDtype == 'AUC'){
        
        cut = seq(0.00000001, 1, length.out = FDcut_number)
        width = diff(cut)
        
        if (datatype == 'abundance') {
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(data, FDdistM, tau, datatype = 'abundance')
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(gamma, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
          alpha = data.frame(alpha, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
          beta = data.frame(beta, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
        }
        
        if (datatype == 'incidence_raw') {
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, datatype = 'incidence_raw')
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(gamma, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
          alpha = data.frame(alpha, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
          beta = data.frame(beta, q) %>% set_colnames(c('Diversity', 'Order.q'))
          
        }
        
      }
      
      C = beta %>% mutate(Diversity = ifelse(Order.q == 1, log(Diversity)/log(N), (Diversity^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      U = beta %>% mutate(Diversity = ifelse(Order.q == 1, log(Diversity)/log(N), (Diversity^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      V = beta %>% mutate(Diversity = (Diversity - 1)/(N - 1))
      S = beta %>% mutate(Diversity = (1/Diversity - 1)/(1/N - 1))
      
      if(nboot > 1){
        
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)
            
            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, 'abundance')
            
            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
            
            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector
            
            if (FDtype == 'tau_value'){
              
              output = FD_by_tau(data_bt, distance_matrix_bt, FDtau, datatype='abundance')
              gamma = output$gamma
              alpha = output$alpha
              beta=gamma/alpha
              
            }
            
            if (FDtype == 'AUC'){
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(data_bt, distance_matrix_bt, tau, datatype = 'abundance')
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              beta_over_tau = gamma_over_tau/alpha_over_tau
              
              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              
              beta = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          if (datatype == 'incidence_raw') {
            
            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])
            
            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), FDdistM, f0_hat, 'incidence_freq')
            
            raw = lapply(1:ncol(p_bt), function(j){
              
              lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))
            
            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)
            
            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
            
            if (FDtype == 'tau_value'){
              
              output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, datatype = 'incidence_raw')
              gamma = output$gamma
              alpha = output$alpha
              beta = gamma/alpha
              
            }
            
            if (FDtype == 'AUC'){
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, datatype = 'incidence_raw')
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              beta_over_tau = gamma_over_tau/alpha_over_tau
              
              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              
              beta = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          beta = gamma/alpha
          
          beta = data.frame(Diversity = beta, Order.q = q)
          
          C = (beta %>% mutate(Diversity = ifelse(Order.q == 1, log(Diversity)/log(N), (Diversity^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Diversity
          U = (beta %>% mutate(Diversity = ifelse(Order.q == 1, log(Diversity)/log(N), (Diversity^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Diversity
          V = (beta %>% mutate(Diversity = (Diversity - 1)/(N - 1)))$Diversity
          S = (beta %>% mutate(Diversity = (1/Diversity - 1)/(1/N - 1)))$Diversity
          
          beta = beta$Diversity
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }) %>% abind(along = 3) %>% apply(1:2, sd)
        
      } else {
        
        se = matrix(0, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    
    se = as.data.frame(se)
    
    if (diversity == "TD") index = "TD"
    if (diversity == "PD" & PDtype == "PD") index = "PD"
    if (diversity == "PD" & PDtype == "meanPD") index = "meanPD"
    if (diversity == "FD" & FDtype == "tau_value") index = "FD_tau"
    if (diversity == "FD" & FDtype == "AUC") index = "FD_AUC"
    
    gamma = gamma %>% mutate(s.e. = se$gamma,
                             LCL = Diversity - tmp * se$gamma,
                             UCL = Diversity + tmp * se$gamma,
                             Region = region_name,
                             diversity = index)
    
    alpha = alpha %>% mutate(s.e. = se$alpha,
                             LCL = Diversity - tmp * se$alpha,
                             UCL = Diversity + tmp * se$alpha,
                             Region = region_name,
                             diversity = index)
    
    beta = beta %>% mutate(  s.e. = se$beta,
                             LCL = Diversity - tmp * se$beta,
                             UCL = Diversity + tmp * se$beta,
                             Region = region_name,
                             diversity = index)
    
    C = C %>% mutate(        s.e. = se$C,
                             LCL = Diversity - tmp * se$C,
                             UCL = Diversity + tmp * se$C,
                             Region = region_name,
                             diversity = index)
    
    
    U = U %>% mutate(        s.e. = se$U,
                             LCL = Diversity - tmp * se$U,
                             UCL = Diversity + tmp * se$U,
                             Region = region_name,
                             diversity = index)
    
    V = V %>% mutate(        s.e. = se$V,
                             LCL = Diversity - tmp * se$V,
                             UCL = Diversity + tmp * se$V,
                             Region = region_name,
                             diversity = index)
    
    S = S %>% mutate(        s.e. = se$S,
                             LCL = Diversity - tmp * se$S,
                             UCL = Diversity + tmp * se$S,
                             Region = region_name,
                             diversity = index)
    
    rbind(gamma, alpha, beta, C, U, V, S) %>% cbind(., type = rep(c('gamma', 'alpha', 'beta', 'C', 'U', 'V', 'S'), each = length(q)))
    
  }
  
  output = lapply(1:length(data_list), function(i) for_each_region(data = data_list[[i]], region_name = region_names[i], N = Ns[i]))
  names(output) = region_names
  
  return(output)
  
}



#' ggplot for beta diversity
#' 
#' \code{ggobsbeta3D}: ggplot for observed beta diversity with order q
#' 
#' @param output the output from obsbeta3D
#' @param type selection of plot type : \code{type = 'B'} for plotting the gamma, alpha, and beta diversity ; \code{type = 'D'} for plotting 4 turnover dissimilarities.
#' @param measurement character indicating the label of y-axis.
#' @param transp a value between 0 and 1 controlling transparency. \code{transp = 0} is completely transparent, default is 0.4.
#' 
#' @return a figure for Beta diversity or dissimilarity diversity.
#' 
#' @examples
#' ## Taxonomic diversity for abundance data
#' data(beetle_abu)
#' output1 = obsbeta3D(data = beetle_abu, diversity = 'TD', datatype = 'abundance', nboot = 20)
#' 
#' ggobsbeta3D(output1, type = 'B')
#' ggobsbeta3D(output1, type = 'D')
#' 
#' 
#' ## Taxonomic diversity for incidence data
#' data(beetle_inc)
#' output2 = obsbeta3D(data = beetle_inc, diversity = 'TD', datatype = 'incidence_raw', nboot = 20, conf = 0.95)
#' 
#' ggobsbeta3D(output2, type = 'B')
#' ggobsbeta3D(output2, type = 'D')
#' 
#' 
#' ## Phylogenetic Hill numbers for abundance data
#' data(beetle_abu)
#' data(beetle_tree)
#' output3 = obsbeta3D(data = beetle_abu, diversity = 'PD', datatype = 'abundance', nboot = 20, conf = 0.95, PDtree = beetle_tree, PDreftime = NULL, PDtype = 'meanPD')
#' 
#' ggobsbeta3D(output3, type = 'B')
#' ggobsbeta3D(output3, type = 'D')
#' 
#' 
#' ## Phylogenetic Hill numbers for incidence data
#' data(beetle_inc)
#' data(beetle_tree)
#' output4 = obsbeta3D(data = beetle_inc, diversity = 'PD', datatype = 'incidence_raw', nboot = 0, conf = 0.95, PDtree = beetle_tree, PDreftime = NULL, PDtype = 'meanPD')
#' 
#' ggobsbeta3D(output4, type = 'B')
#' ggobsbeta3D(output4, type = 'D')
#' 
#' 
#' ## Functional diversity for abundance data under single threshold
#' data(beetle_abu)
#' data(beetle_distM)
#' output5 = obsbeta3D(data = beetle_abu, diversity = 'FD', datatype = 'abundance', nboot = 20, conf = 0.95, FDdistM = beetle_distM, FDtype = 'tau_value', FDtau = NULL)
#' 
#' ggobsbeta3D(output5, type = 'B')
#' ggobsbeta3D(output5, type = 'D')
#' 
#' 
#' ## Functional diversity for incidence data under single threshold
#' data(beetle_inc)
#' data(beetle_distM)
#' output6 = obsbeta3D(data = beetle_inc, diversity = 'FD', datatype = 'incidence_raw', nboot = 20, conf = 0.95, FDdistM = beetle_distM, FDtype = 'tau_value', FDtau = NULL)
#'
#' ggobsbeta3D(output6, type = 'B')
#' ggobsbeta3D(output6, type = 'D')
#' 
#' 
#' ## Functional diversity for abundance data with thresholds integrating from 0 to 1
#' data(beetle_abu)
#' data(beetle_distM)
#' output7 = obsbeta3D(data = beetle_abu, diversity = 'FD', datatype = 'abundance', nboot = 10, conf = 0.95, FDdistM = beetle_distM, FDtype = 'AUC', FDcut_number = 30)
#' 
#' ggobsbeta3D(output7, type = 'B')
#' ggobsbeta3D(output7, type = 'D')
#' 
#' 
#' ## Functional diversity for incidence data with thresholds integrating from 0 to 1
#' data(beetle_inc)
#' data(beetle_distM)
#' output8 = obsbeta3D(data = beetle_inc, diversity = 'FD', datatype = 'incidence_raw', nboot = 10, conf = 0.95, FDdistM = beetle_distM, FDtype = 'AUC', FDcut_number = 30)
#' 
#' ggobsbeta3D(output8, type = 'B')
#' ggobsbeta3D(output8, type = 'D')
#' 
#' @export
ggobsbeta3D = function(output, type = c('B', 'D')){
  
  if (type == 'B'){
    
    gamma = lapply(output, function(y) y %>% filter(type == 'gamma')) %>% do.call(rbind,.) %>% mutate(div_type = "Gamma") %>% as_tibble()
    alpha = lapply(output, function(y) y %>% filter(type == 'alpha')) %>% do.call(rbind,.) %>% mutate(div_type = "Alpha") %>% as_tibble()
    beta =  lapply(output, function(y) y %>% filter(type == 'beta'))  %>% do.call(rbind,.) %>% mutate(div_type = "Beta")  %>% as_tibble()
    
    df = rbind(gamma, alpha, beta)
    df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))
    
  }
  
  if (type == 'D'){
    
    C = lapply(output, function(y) y %>% filter(type == 'C')) %>% do.call(rbind,.) %>% mutate(div_type = "1-CqN") %>% as_tibble()
    U = lapply(output, function(y) y %>% filter(type == 'U')) %>% do.call(rbind,.) %>% mutate(div_type = "1-UqN") %>% as_tibble()
    V = lapply(output, function(y) y %>% filter(type == 'V')) %>% do.call(rbind,.) %>% mutate(div_type = "1-VqN") %>% as_tibble()
    S = lapply(output, function(y) y %>% filter(type == 'S')) %>% do.call(rbind,.) %>% mutate(div_type = "1-SqN") %>% as_tibble()
    
    df = rbind(C, U, V, S)
    df$div_type <- factor(df$div_type, levels = c("1-CqN", "1-UqN", "1-VqN", "1-SqN"))
    
  }
  
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  
  if (unique(output[[1]]$diversity) == 'TD'    ) { ylab = "Taxonomic diversity" }
  if (unique(output[[1]]$diversity) == 'PD'    ) { ylab = "Phylogenetic diversity" }
  if (unique(output[[1]]$diversity) == 'meanPD') { ylab = "Phylogenetic Hill number" }
  if (unique(output[[1]]$diversity) == 'FD_tau') { ylab = "Functional diversity (given tau)" }
  if (unique(output[[1]]$diversity) == 'FD_AUC') { ylab = "Functional diversity (AUC)" }
  
  ggplot(data = df, aes(x = Order.q, y = Diversity, col = Region)) +
    geom_line(aes(x = Order.q, y = Diversity, col = Region), size=1.2) + 
    geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Region, col = NULL), alpha = 0.4) + 
    scale_colour_manual(values = cbPalette) + 
    scale_fill_manual(values = cbPalette) + 
    facet_grid(div_type ~ ., scales = 'free_y') +
    theme_bw() + 
    theme(legend.position = "bottom", legend.title = element_blank()) +
    labs(x = 'Order q', y = ylab)
}


bootstrap_population_multiple_assemblage = function(data, data_gamma, datatype){
  
  if (datatype == 'abundance'){
    
    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()
    
    output = apply(data, 2, function(x){
      
      p_i_hat = iNEXT.3D:::EstiBootComm.Ind(Spec = x)
      
      if(length(p_i_hat) != length(x)){
        
        p_i_hat_unobs = p_i_hat[(length(x)+1):length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat==0)
        
        chosen = sample(x = candidate, size = min(length(p_i_hat_unobs), length(candidate)), replace = F)
        p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)
        
        p_i_hat
        
      } else {
        
        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat
        
      }
    })
    
  }
  
  if (datatype == 'incidence'){
    
    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling
    
    output = apply(data, 2, function(x){
      
      pi_i_hat = iNEXT.3D:::EstiBootComm.Sam(Spec = x)
      
      if(length(pi_i_hat) != (length(x) - 1)){
        
        pi_i_hat_unobs = pi_i_hat[length(x):length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:(length(x)-1)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat == 0)
        chosen = sample(x = candidate, size = min(length(pi_i_hat_unobs), length(candidate)), replace = F)
        pi_i_hat[chosen] = pi_i_hat_unobs
        
        pi_i_hat
        
      } else {
        
        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat
        
      }
    })
    
  }
  
  return(output)
  
}

Bootstrap_distance_matrix = function(data, distance_matrix, f0.hat, datatype){
  
  if (datatype == "incidence_freq") {
    n = data[1]
    X = data[-1]
    u = sum(data)
  } else if (datatype == "abundance") {
    n = sum(data)
    X = data
  }
  
  # n = sum(data)
  distance = as.matrix(distance_matrix)
  dij = distance
  # X = data
  
  F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
  F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
  
  if (datatype == "abundance") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  } else if (datatype == "incidence_freq") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
  }
  
  if (f0.hat == 0) {
    d = dij
  } else if (f0.hat == 1) {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  } else {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    fo.num = (f0.hat * (f0.hat-1) )/2
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)/fo.num
    d00 <- pmax(d00, t(d00))###signmatrix
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  }
  
  return(d)
  
}

FD.m.est_0 = function (ai_vi, m, q, nT) {
  EFD = function(m, qs, obs, asy, beta, av) {
    m = m - nT
    out <- sapply(1:length(qs), function(i) {
      if (qs[i] != 2) {
        obs[i] + (asy[i] - obs[i]) * (1 - (1 - beta[i])^m)
      }
      else if (qs[i] == 2) {
        V_bar^2/sum((av[, 2]) * ((1/(nT + m)) * (av[, 1]/nT) + ((nT + m - 1)/(nT + m)) * (av[, 1] * (av[, 1] - 1)/(nT * (nT - 1)))))
      }
    })
    return(out)
  }
  V_bar <- sum(ai_vi$ai[, 1] * ai_vi$vi[, 1])/nT
  asy <- iNEXT.3D:::FD_est(ai_vi, q, nT)$est
  obs <- iNEXT.3D:::FD_mle(ai_vi, q)
  out <- sapply(1:ncol(ai_vi$ai), function(i) {
    ai <- ai_vi$ai[, i]
    ai[ai < 1] <- 1
    av = cbind(ai = round(ai), vi = ai_vi$vi[, i])
    RFD_m = iNEXT.3D:::RFD(av, nT, nT - 1, q, V_bar)
    beta <- rep(0, length(q))
    asy_i <- asy[, i]
    obs_i <- obs[, i]
    asy_i <- sapply(1:length(q), function(j) {
      max(asy_i[j], obs_i[j], RFD_m[j])
    })
    
    obs_i <- sapply(1:length(q), function(j) {
      max(RFD_m[j], obs_i[j])
    })
    
    beta0plus <- which(asy_i != obs_i)
    beta[beta0plus] <- (obs_i[beta0plus] - RFD_m[beta0plus])/(asy_i[beta0plus] - RFD_m[beta0plus])
    
    if (sum(m < nT) != 0) {
      int.m = sort(unique(c(floor(m[m < nT]), ceiling(m[m < nT]))))
      mRFD = rbind(int.m, sapply(int.m, function(k) iNEXT.3D:::RFD(av, nT, k, q, V_bar)))
    }
    
    sapply(m, function(mm) {
      if (mm < nT) {
        if (mm == round(mm)) {
          mRFD[-1, mRFD[1, ] == mm]
        } else {
          (ceiling(mm) - mm) * mRFD[-1, mRFD[1, ] == floor(mm)] + (mm - floor(mm)) * mRFD[-1, mRFD[1, ] == ceiling(mm)]
        }
      }
      else if (mm == nT) {
        obs_i
      }
      else if (mm == Inf) {
        asy_i
      }
      else {
        EFD(m = mm, qs = q, obs = obs_i, asy = asy_i, beta = beta, av = av)
      }
    }) %>% t %>% as.numeric
  })
  matrix(out, ncol = ncol(ai_vi$ai))
}
