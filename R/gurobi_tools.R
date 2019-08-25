# library(gurobi)
# library(sybil)
# library(dplyr)

insert_into_list <- function(list, insert, idx){
  segment_1 <- list[1:idx]
  segment_2 <- list[(idx+1):length(list)]
  
  new_list <- c(segment_1, insert, segment_2)
  return(new_list)
}

convert_sybil_to_lp <- function(sybil, output = 'temp.lp'){
  prob <- sysBiolAlg(sybil, useNames=TRUE)
  writeProb(problem(prob), fname = output)
}

#' description
#' @param model filename of RData file for sybil model; defaults to 'mutans_model.RData'
#' @return GRB model
#' @seealso 
#' @export
#' @examples
#' 
convert_sybil_to_gurobi <- function(sybil, output = 'temp.lp'){
  convert_sybil_to_lp(sybil, output = output)
  model <- gurobi_read(output)
  return(model)
}

# MODELS


#' description
#' @return GRB model
#' @seealso 
#' @export
#' @examples
#' 
GRB_ecoli_model <- function(){
  data(Ec_core);
  model=Ec_core;
  model <- changeBounds(model, 11, lb = 0)
  # model <- changeBounds(model, 13, lb = 0, ub = 0)
  model <- rmReact(model = model, react = 13)
  for (i in findExchReact(model)@react_pos){
    model <- changeBounds(model, i, lb = -1000, ub = 1000)
    # if (model@lowbnd[i] == 0){
    #   model <- changeBounds(model, i, lb = -1000)
    # }
  }
  
  ecoli <- convert_sybil_to_gurobi(model)
  # ecoli$show_output(FALSE)
  return(ecoli)
}

#' description
#' @return GRB falcon model
#' @seealso 
#' @export
#' @examples
#' 
GRB_ecoli_falcon_model <- function(){
  sybil_ecoli <- get_ecoli_model()
  ecoli_falcon_model <- GRB_generate_falcon_model(sybil_ecoli)
  return(ecoli_falcon_model)
}

#' Gurobi S. cerevisiae model
#' @param model filename of RData file for sybil model; defaults to 'yeast_model.RData'
#' @return GRB model
#' @seealso 
#' @export
#' @examples
#' 
GRB_yeast_model <- function(model = 'yeast_model.RData'){
  #maranas_model_exch_add_biom_rm/maranas_model_lipid_exch
  load(model)
  
  yeast <- convert_sybil_to_gurobi(yeast_model)
  # yeast$show_output(FALSE)
  
  return(yeast)
}

#' Gurobi S. cerevisiae falcon model
#' @param model filename of RData file for sybil model; defaults to 'yeast_model.RData'
#' @return GRB falcon model
#' @seealso 
#' @export
#' @examples
#' 
GRB_yeast_falcon_model <- function(model = 'yeast_model.RData'){
  
  load(model)
  
  sybil_yeast <- yeast_model
  yeast_falcon_model <- GRB_generate_falcon_model(sybil_yeast)
  return(yeast_falcon_model)
}

#' Gurobi S. mutans model
#' @param model filename of RData file for sybil model; defaults to 'mutans_model.RData'
#' @return GRB model
#' @seealso 
#' @export
#' @examples
#' 
GRB_mutans_model <- function(model = 'mutans_model.RData'){
  
  load(model)
  
  mutans <- convert_sybil_to_gurobi(mutans_model)
  # mutans$show_output(FALSE)
  
  return(mutans)
}

#' Gurobi S. mutans falcon model
#' @param model filename of RData file for sybil model; defaults to 'mutans_model.RData'
#' @return GRB falcon model
#' @seealso 
#' @export
#' @examples
#' 
GRB_mutans_falcon_model <- function(model = 'mutans_model.RData'){
  load(model)
  
  sybil_mutans <- mutans
  mutans_falcon_model <- GRB_generate_falcon_model(sybil_mutans)
  return(mutans_falcon_model)
}

#' Gurobi P. aeruginosa model
#' @param model filename of RData file for sybil model; defaults to 'pao_model.RData'
#' @return GRB model
#' @seealso 
#' @export
#' @examples
#' 
GRB_pao_model <- function(model = 'pao_model.RData'){
  load(model)
  
  pao <- convert_sybil_to_gurobi(pao_model)
  # pao$show_output(FALSE)
  
  return(pao)
}

#' Gurobi P. aeruginosa falcon model
#' @param model filename of RData file for sybil model; defaults to 'pao_model.RData'
#' @return GRB falcon model
#' @seealso 
#' @export
#' @examples
#' 
GRB_pao_falcon_model <- function(model = 'pao_model.RData'){
  load(model)
  
  pao_falcon <- GRB_generate_falcon_model(pao_model)
  # pao_falcon$show_output(FALSE)
  
  return(pao_falcon)
}


# FUNCTIONS

#' Generate Gurobi falcon model. First generates model in sybil then converts the model to Gurobi in order to add additional
#' constraints on flux
#' @param sybil_model base model to build falcon model from
#' @param falcon_model boolean value indicating whether or not the base model has enzyme activity already. If TRUE, then the function will simply add contstraints to the falcon model
#' @param r0_gene_set list of genes passed in to the falcon-model building method (generate_falcon_model)
#' @param r0_rxn_set_list list of reaction sets passed in to the falcon-model building method (generate_falcon_model)
#' @return Gurobi model based on metabolic model with enzymatic acctivity
#' @seealso 
#' @export
#' @examples
#' 
GRB_generate_falcon_model <- function(sybil_model, falcon_model = FALSE, r0_gene_set = c(), r0_rxn_set_list = c()){
  
  sybil_falcon_model <- sybil_model
  
  if (!falcon_model){ # if user has not passed in a falcon model already
    sybil_falcon_model <- generate_falcon_model(sybil_model, r0_gene_set, r0_rxn_set_list)
  }
  
  ## ADD NECESSARY CONSTRAINTS TO MODEL
  vars <- sybil_falcon_model@react_id #grb_falcon_model$get_names()$VarName
  
  split_fwd_rxns <- vars[grep('fwd', vars)]
  split_rev_rxns <- vars[grep('rev', vars)]
  
  split_rxns <- sapply(split_fwd_rxns, function(x) strsplit(x, ' ')[[1]][1])
  split_rxns <- unique(split_rxns)
  
  if (length(split_fwd_rxns) != length(split_rev_rxns)){
    print('fed rev matchup error')
    return()
  }
  
  Binaries <- c('', 'Binaries')
  
  for (rxn in split_rxns){
    I <- paste('I', rxn, sep = '_')
    Binaries <- c(Binaries, I)
  }
  
  new_constr <- c()
  
  # add bounds on split conversion reactions (gene -> activity_[rxn])
  for (i in 1:length(split_fwd_rxns)){
    fwd <- split_fwd_rxns[i]
    fwd <- gsub(pattern = ' ', replacement = '_', fwd)
    rev <- split_rev_rxns[i]
    rev <- gsub(pattern = ' ', replacement = '_', rev)
    
    a_rxn <- strsplit(fwd, ' ')[[1]][1]
    I <- paste('I', a_rxn, sep = '_')
    
    fwd_name <- paste(fwd, 'I', sep = '_')
    rev_name <- paste(rev, 'I', sep = '_')
    
    fwd_line <- paste(' ', fwd_name, ": - ", fwd, ' + 1000 ', I, ' >= 0', sep = '')
    rev_line <- paste(' ',rev_name, ": - ", rev, ' - 1000 ', I, ' >= -1000', sep = '')
    
    new_constr <- c(new_constr, fwd_line, rev_line)
  }
  
  # write and read model
  convert_sybil_to_lp(sybil_falcon_model)
  lp_lines <- rewrite_lp(model = sybil_falcon_model)
  
  constr_idx <- grep('Bounds', lp_lines)
  lp_lines <- insert_into_list(lp_lines, new_constr, constr_idx-2)
  end_idx <- grep('End', lp_lines)
  lp_lines <- lp_lines <- insert_into_list(lp_lines, Binaries, (end_idx-2))
  write(lp_lines, file = 'temp.lp')
  model <- gurobi_read('temp.lp')
  
  return(model)
}

#' Solve Gurobi model
#' @param model Gurobi model
#' @return solution to optimization problem
#' @seealso 
#' @export
#' @examples
#' 
solve_model <- function(model, obj_idx, sense = 'max'){
  model$obj <- rep(0, length(model$obj))
  model$obj[obj_idx] <- 1
  model$modelsense <- sense
  sol <- gurobi(model, list(OutputFlag = 0))
  return(sol)
}

#' determine flux coupling between reactions in a Gurobi model
#' @param model Gurobi model
#' @param min_fva minimum correlation value to test for coupling
#' @param fix_frac const value used in fixing flux at non-zero value 
#' @param fix_tol_frac maximum allowed error when determining whether flux for a reaction is fixed
#' @param bnd_tol allowed error in comparing max & min flux
#' @param stored_obs maximum number of flux instancecs to store
#' @param cor_iter number of iterations after which correlation is considered in checking coupling
#' @param reaction_indexes indexes for reactions to check; check all if empty
#' @param compare_mtx boolean to indicate whether or not to couple known sets together if any individual reactions are coupled; FALSE by default
#' @param known_set_mtx boolean matrix of known coupling between reactions
#' @param directional boolean to indicate whether to check for directional coupling between reactions
#' @return boolean matrix indicating whether reactions are coupled to each other (reactions are assigned to rows and columns)
#' @seealso 
#' @export
#' @examples
#' 
flux_coupling_raptor <- function(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01,
                                 bnd_tol = 0.1, stored_obs = 4000, cor_iter = 5, cor_check = TRUE,
                                 reaction_indexes = c(), compare_mtx = FALSE, 
                                 known_set_mtx = Matrix(data = FALSE, nrow = 1, ncol = 1, sparse = TRUE),
                                 directional = FALSE) {
  
  # min_fva_cor is minimum correlation between fluxes
  # bnd_tol is allowed error in comparing max & min flux
  # fix_frac is const value used in fixing flux at non-zero value
  # fix_tol_frac is error allowed in determining whether flux is fixed
  # stored_obs is # not flux values to be stored
  # cor_iter is number of iterations after which correlation is considered in checking coupling
  
  if (!cor_check){
    stored_obs = 1
  }
  
  vars <- model$varnames
  n <- length(vars)
  
  if (is.null(known_set_mtx)){
    compare_mtx <- FALSE
  }
  else {
    if ((nrow(known_set_mtx) < n) | (ncol(known_set_mtx) < n)){compare_mtx <- FALSE}
  }
  
  # if empty set, then assume all reactions are to be inspected
  if (length(reaction_indexes) == 0){
    reaction_indexes <- c(1:n)
  }
  
  prev_obj <- model$obj
  model$obj <- rep(0,n) # clear the objective
  prev_sense <- model$modelsense
  
  original_ub <- model$ub
  original_lb <- model$lb
  
  global_max <- rep(0, n)
  global_min <- rep(0, n)
  
  coupled <- Matrix(FALSE, nrow=n, ncol=n, dimnames=list(vars, vars), sparse = TRUE)
  
  blocked <- rep(FALSE, n)
  active <- rep(TRUE, n)
  
  not_fixed <- function(x,y) { # check for variability in flux, return TRUE or FALSE
    !is.infinite(x) && !is.infinite(y) && abs(x - y) > fix_tol_frac*max(abs(x), abs(y))
  }
  
  rxn_fix <- function(max_, min_){
    
    avg <- min_ + fix_frac*(max_ - min_)
    ct = 0
    while (near(avg, 0) & ct < 5){
      print(c(avg, max_, min_))
      avg <- avg  + fix_frac*(max_ - min_)
      ct = ct + 1
    }
    
    return(avg)
  }
  
  correlation_check <- function(flux, i, j){ # return true if correlation is high, or no correlation --> continue comparing flux
    # if false, skip in depth comparison
    
    n_entries <- length(which(!near(flux[,i], 0, tol = bnd_tol) | !near(flux[,j], 0, tol = bnd_tol)))
    if (n_entries > cor_iter){
      C <- cor(flux[,i], flux[,j])
      
      if ((is.na(C)) | (abs(C) < min_fva_cor)){
        return(FALSE)
      }
    }
    return(TRUE)
  }
  
  flux <- matrix(c(0), nrow = stored_obs, ncol = n)
  lp_calls <- 0
  
  update_flux <- function(flux_, idx, sol){
    if (stored_obs > 0){
      flux_[idx,] <- sol
    }
    return(flux_)
  }
  
  # if j is coupled to i, couple the rxns known to be coupled to j to i as well
  known_set_coupling <- function(i, j, coupled, active){
    set <- which(known_set_mtx[j,])
    
    coupled[i, set] <- TRUE
    coupled[set, set] <- TRUE
    active[set] <- FALSE
    
    list(coupled = coupled, active = active)
  }
  
  check_directional_coupling <- function(i, j, fixed = TRUE){
    
    info <- list(lp_calls = 0,
                 coupled = FALSE)
    
    prev_ub_i <- model$ub[i]
    prev_lb_i <- model$lb[i]
    
    if (!fixed){
      fixed_val <- 0.01
      if (!near(global_max[i], 0, tol = bnd_tol) | !near(global_min[i], 0, tol = bnd_tol)){
        fixed_val <- rxn_fix(global_max[i], global_min[i])
      }
      else {
        if (!near(model$ub[i], 0)){ #model$getattr("UB")[vars[i]] > tol
          
          model$obj[i] <- 1
          model$modelsense <- 'max'
          sol <- gurobi(model, list(OutputFlag = 0))
          info$lp_calls <- info$lp_calls + 1
          global_max <- pmax(global_max, sol$x)
          global_min <- pmin(global_min, sol$x)
          #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
          # flux[lp_calls%%stored_obs,] <- sol$x
          
          if (!near(global_max[i], 0) | !near(global_min[i], 0)){
            fixed_val <- rxn_fix(global_max[i], global_min[i])
          }
          
          model$obj <- rep(0, n)
        }
        if (!near(model$lb[i], 0)){ #model$getattr("LB")[vars[i]] < (-1*tol_)
          
          model$obj[i] <- 1
          model$modelsense <- 'min'
          sol <- gurobi(model, list(OutputFlag = 0))
          info$lp_calls <- info$lp_calls + 1
          #print(is.null(flux))
          global_max <- pmax(global_max, sol$x)
          global_min <- pmin(global_min, sol$x)
          #flux <- update_flux(flux, lp_calls%%stored_obs, sol$X)
          # flux[lp_calls%%stored_obs,] <- sol$x
          
          if (!near(global_max[i], 0) | !near(global_min[i], 0)){
            fixed_val <- rxn_fix(global_max[i], global_min[i])
          }
          
          model$obj <- rep(0, n)
        }
        
        model$ub[i] <- fixed_val
        model$lb[i] <- fixed_val
      }
    }
    
    prev_ub_j <- model$ub[j]
    prev_lb_j <- model$lb[j]
    
    # zero out target bounds
    model$ub[j] <- 0
    model$lb[j] <- 0
    
    # check feasibility
    sol <- gurobi(model, list(OutputFlag = 0))
    
    info$lp_calls <- info$lp_calls + 1
    
    feasible <- !(length(sol$x) == 0)
    info$coupled <- feasible
    
    model$ub[j] <- prev_ub_j
    model$lb[j] <- prev_lb_j
    
    model$ub[i] <- prev_ub_i
    model$lb[i] <- prev_lb_i
    
    
    return(info)
  }
  
  
  for (idx in 1:(length(reaction_indexes))) { # (i in 1:(n-1))
    # iterate over passed in idxs instead (idx in 1:length(reaction_indexes)); i <-  reaction_indexes[idx]
    i <-  reaction_indexes[idx]
    
    if (!active[i] | blocked[i]) next
    sub_max <- rep(-Inf, n)
    sub_min <- rep(Inf, n)
    prev_ub <- model$ub[i]
    prev_lb <- model$lb[i]
    
    if (!near(global_max[i], 0, tol = bnd_tol) | !near(global_min[i], 0, tol = bnd_tol)){
      fixed_val <- rxn_fix(global_max[i], global_min[i])
    }
    else {
      if (!near(model$ub[i], 0)){
        
        model$obj[i] <- 1
        model$modelsense <- 'max'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        if (!near(global_max[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }
        
        model$obj <- rep(0, n)
      }
      if (!near(model$lb[i], 0)){
        
        model$obj[i] <- 1
        model$modelsense <- 'min'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        flux[lp_calls%%stored_obs,] <- sol$x
        
        if (!near(global_max[i], 0) | !near(global_min[i], 0)){
          fixed_val <- rxn_fix(global_max[i], global_min[i])
        }
        
        model$obj <- rep(0, n)
      }
      if (near(global_max[i], 0) & near(global_min[i], 0)){ #(abs(global_max[i]) < tol) & (abs(global_min[i]) < tol)
        blocked[i] <- TRUE
        active[i] <- FALSE
        next
      }
    }
    
    # set new bounds for selected rxn (temporarily)
    model$ub[i] <- fixed_val
    model$lb[i] <- fixed_val
    
    # couple reaction to itself if not blocked
    if (!blocked[i]){
      coupled[i,i] <- TRUE
      active[i] <- FALSE
    }
    
    if (idx == length(reaction_indexes)){break}
    
    for (idx2 in (idx+1):length(reaction_indexes)) { # (j in (i+1):n)
      # also keep this in passed in idxs (idx2 in (idx+1):length(reaction_indexes)); j <-  reaction_indexes[idx2]
      j <-  reaction_indexes[idx2]
      # check for fixed or blocked
      if (!active[j] | blocked[j]){next}
      if (not_fixed(sub_max[j], sub_min[j])){next}

      # check for uncoupled via correlation
      if (cor_check){
        if (!correlation_check(flux, i, j)){next}
      }
      
      skip <- FALSE
      
      max <- 0
      min <- 0
      
      model$obj[j] <- 1
      if (!skip) {
        model$modelsense <- 'max'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        sub_max <- pmax(sub_max, sol$x)
        sub_min <- pmin(sub_min, sol$x)
        
        flux[lp_calls%%stored_obs,] <- sol$x
        
        max <- sol$x[j]
        
        skip <- not_fixed(sub_max[j], sub_min[j])
      }
      
      if (!skip) {
        model$modelsense <- 'min'
        sol <- gurobi(model, list(OutputFlag = 0))
        lp_calls <- lp_calls + 1
        global_max <- pmax(global_max, sol$x)
        global_min <- pmin(global_min, sol$x)
        sub_max <- pmax(sub_max, sol$x)
        sub_min <- pmin(sub_min, sol$x)
        
        flux[lp_calls%%stored_obs,] <- sol$x
        
        max <- sol$x[j]
        
        skip <- not_fixed(sub_max[j], sub_min[j])
      }
      
      if (near(max, 0) & near(min, 0)){skip = TRUE}
      
      if (!skip) { # finally label as coupled
        coupled[i,j] <- TRUE
        coupled[j,j] <- TRUE # make sure the reaction doesn't look blocked
        active[j] <- FALSE
        
        if (compare_mtx){
          output <- known_set_coupling(i, j, coupled, active)
          coupled <- output$coupled
          active <- output$active
        }
      }
      
      if (skip & directional){
        dir_i_j <- check_directional_coupling(i, j, fixed = TRUE)
        lp_calls <- lp_calls + dir_i_j$lp_calls
        
        if (dir_i_j$coupled){
          coupled[i,j] <- TRUE
          next
        }
        
        dir_j_i <- check_directional_coupling(j, i, fixed = FALSE)
        lp_calls <- lp_calls + dir_j_i$lp_calls
        
        if (dir_j_i$coupled){
          coupled[j,i] <- TRUE
        }
        
        if (dir_i_j$coupled | dir_j_i$coupled){
          coupled[i,j] <- TRUE
          coupled[j,i] <- TRUE
          active[j] <- FALSE
        }
      }
      
      
      model$obj <- rep(0, n)
    }
    
    
    # unfix i
    model$ub <- original_ub
    model$lb <- original_lb
  }
  
  model$obj <- prev_obj
  model$modelsense <- prev_sense
  
  print(lp_calls)
  
  list(
    coupled = coupled,
    lp_calls = lp_calls
  )
}

#' get the position of a reaction along its axis on the S matrix
#' @param model Gurobi model
#' @param rxn reaction id to get the index of
#' @return integer indicating the index of the reaction along the S matrix
#' @seealso 
#' @export
#' @examples
#' 
GRB_get_rxn_idx <- function(model, rxn){
  vars <- model$varnames
  return(return(which(vars == rxn)))
}

GRB_generate_pair_list <- function(model_og){
  model <- model_og
  return(return_couples(flux_coupling_raptor(model)$coupled))
}

GRB_generate_set_list <- function(model_og, reaction_indexes = c()){
  model <- model_og
  return(get_list_of_sets(return_couples(flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled)))
}

#' wrapper function for flux coupling
#' @param i integer indicating index of the reaction to suppress
#' @param vars vector of reaction ids
#' @param model_og Gurobi model
#' @param reaction_indexes integers indicating reactions to run flux coupling over; defaults to 1:length(vars)
#' @param compare_mtx boolean to indicate whether or not to couple known sets together if any individual reactions are coupled; FALSE by default
#' @param init_coupling_mtx boolean matrix of known coupling between reactions; empty by default
#' @param file_output filename to output coupling vector to
#' @return list of integers which indicate positions of TRUE values in couling matrix
#' @seealso 
#' @examples
#' 
GRB_flux_coupling_raptor_wrapper <- function(i, vars, model_og, reaction_indexes = 1:length(vars), 
                                             compare_mtx = FALSE, init_coupling_mtx = c(), file_output = NULL){
  print(paste('suppression index:', i))
  
  model <- model_og
  
  # block i
  model$setattr("UB", setNames(0, vars[i]))
  model$setattr("LB", setNames(0, vars[i]))
  
  coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes, compare_mtx = compare_mtx, 
                                       known_set_mtx = init_coupling_mtx, stored_obs = 4000)$coupled
  
  coupling_idxs <- which(coupling_mtx)
  if (!is.null(file_output)){
    write(paste(c(i,coupling_idxs), collapse = ','),file = file_output,append=TRUE)
  }
  return(coupling_idxs)
}

#' Function for PACT method
#' @param model_og model upon which to run PACT
#' @param suppression_idxs list of integers indicating indexes of reactions to iteratively block and calculate coupling for. DEFAULT is -1, which indicates that all reactions should be suppressed
#' @param reaction_indexes list of integers indicating indexes of reaction to check coupling for. DEFAULT is empty list, which indicates that all reactions in the model should be considered
#' @param compare_known_init_sets boolean indicating whether or not to calculate R0 sets in order to run further optimizations during flux coupling and PACT. DEFAULT: FALSE
#' @param optimize_suppr boolean indicating whether or not to optimize the reactions to suppress based on R0 sets. DEFAULT: FALSE
#' @param optimize_rxns boolean indicating whether or not to optimize the reactions couple. DEFAULT: FALSE
#' @param cores integer indicating the number of cores to run on, if intending to parallelize
#' @param avoid_idxs list of integers indicating indexes to specifically avoid suppressing during PACT method
#' @param file_output filename to output coupling vector to
#' @return a coupling list: list of list of integers which each indicate the couplings induced by each suppression. The index position of each list indicates the index of the reaction which was suppressed and the integers in eacch list indicate the positions of TRUE values in the coupling array
#' @seealso 
#' @export
#' @examples
#' 
GRB_generate_set_lists_cluster <- function(model_og, suppression_idxs = -1, reaction_indexes = c(),
                                           compare_known_init_sets = FALSE, optimize_suppr = FALSE, 
                                           optimize_rxns = FALSE, cores = 1, avoid_idxs = c(), file_output = NULL){
  
  vars <- model_og$varnames
  n <- length(vars)
  
  if (suppression_idxs[1] == -1){
    if (length(reaction_indexes) > 0){
      suppression_idxs = reaction_indexes
    }
    else {
      suppression_idxs = 1:n
    }
  }
  
  # dim: rxns_row, rxns_col, deletions
  model <- model_og
  
  #init_coupling_mtx <- c()
  if (compare_known_init_sets){
    init_coupling_mtx <- flux_coupling_raptor(model, reaction_indexes = reaction_indexes)$coupled
    init_coupling_mtx <- fill_coupling_matrix(init_coupling_mtx)
  }
  suppr_vector <- Matrix(data = FALSE, nrow = 1, ncol = n, sparse = TRUE)
  suppr_vector[suppression_idxs] <- TRUE
  if (compare_known_init_sets & optimize_suppr){
    i <- 1
    while (i <= n){
      if (suppr_vector[i]){ # if tagged to be suppressed
        set_idx <- which(init_coupling_mtx[,i])[1] # which is first reaction (row) i is coupled to
        if (!is.na(set_idx)){
          rxn_idxs <- which(init_coupling_mtx[set_idx,]) # other reactions in set
          # only suppress first reaction in set since, theoretically, suppressing any should have the same effect
          suppr_vector[rxn_idxs] <- FALSE
          suppr_vector[rxn_idxs[1]] <- TRUE
        }
        else {
          suppr_vector[i] <- FALSE
        }
        
      }
      i <- i+1
    }
    
    if (optimize_rxns){
      reaction_indexes <- which(suppr_vector)
    }
  }
  
  if (length(avoid_idxs) > 0){
    suppr_vector[avoid_idxs] <- FALSE
  }
  
  print(paste("# of suppressions:", length(which(suppr_vector)), sep = " "))
  coupling <- mclapply(which(suppr_vector), 
                       function(x) GRB_flux_coupling_raptor_wrapper(x, vars, 
                                                                    model_og, reaction_indexes = reaction_indexes, 
                                                                    compare_mtx = compare_known_init_sets, 
                                                                    init_coupling_mtx = init_coupling_mtx, file_output = file_output), 
                       mc.cores = cores)
  
  return(coupling)
}

#' flux-balance optimization function for Gurobi
#' @param model_og Gurobi model
#' @param obj integer indicating the index of the objective reaction
#' @param suppress list of integers indicating which reactions to suppress before optimizing
#' @param max boolean indicating whether to maximize or minimize the objective. TRUE indicates maximization, FALSE indicates minimization. DEFAULT: TRUE
#' @return
#' @seealso 
#' @export
#' @examples
#' 
GRB_maximize <- function(model_og, obj, suppress = c(), max = TRUE){ # suppress is characters
  model <- model_og
  
  vars <- model$varnames
  n <- length(vars)
  
  # clear obj
  model$obj <- rep(0.0, times = n) 
  
  # set suppressions
  if (length(suppress) > 0){
    suppr_idxs <- which(vars %in% suppress)
    model$ub[suppr_idxs] <- 0
    model$lb[suppr_idxs] <- 0
  }
  
  # set obj
  sense <- 'max'
  if (!max){
    sense <- 'min'
  }
  sol <- solve_model(model = model, obj_idx = obj, sense = sense)
  obj_max <- sol$objval
  return(obj_max)
}
