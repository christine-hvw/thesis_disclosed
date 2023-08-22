
# AUXILIARY FUNCTIONS
# ____________________

# Packages
# library(muxViz)
# library(Matrix)
# library(igraph)
# library(tidymodels)
# library(doParallel)
# library(vip)
# 
# # Test data and parameters
# layers <- 4
# obs <- 100
# 
# node_tensor <- lapply(1:layers, function(x) {
#   graph <- sample_gnp(n = obs, p = .5, directed = FALSE, loops = FALSE)
#   as_adjacency_matrix(graph)
#   }
# )
# 
# layer_tensor <- Matrix(1, layers, layers, sparse = TRUE)
# diag(layer_tensor) <- 0
# 
# M <- BuildSupraAdjacencyMatrixFromEdgeColoredMatrices(node_tensor, layer_tensor, layers, obs)
# 
# Taus <- c(0.2, 0.3, 0.5, 0.1)


# Adaption of muxViz centrality functions ---------------------------------

## Degree centrality ------------------------------------------------------


GetMultiDegree_edit <- function (SupraAdjacencyMatrix, Layers, Nodes, isDirected) {
  
  NodesTensor <- SupraAdjacencyToNodesTensor(binarizeMatrix_edit(SupraAdjacencyMatrix), 
                                             Layers, Nodes)
  
  AggrMatrix <- GetAggregateMatrix(NodesTensor, Layers, Nodes)
  
  MultiInDegreeVector <- Matrix::t(sumR(AggrMatrix, 1))
  
  MultiOutDegreeVector <- sumR(AggrMatrix, 2)
  
  if (!isDirected) {
    MultiDegreeVector <- (MultiInDegreeVector + MultiOutDegreeVector)/2
  }
  else {
    MultiDegreeVector <- MultiInDegreeVector + MultiOutDegreeVector
  }
  return(MultiDegreeVector)
}

# Also edit function that also counts inter-layer edges

GetMultiDegreeSum_edit <- function(SupraAdjacencyMatrix, Layers, Nodes, NodesTensor, isDirected) {
  
  SupraInDegree <- Matrix::t(sumR(binarizeMatrix_edit(SupraAdjacencyMatrix), 1))
  
  MultiInDegreeVector <- sumR(reshapeR(SupraInDegree, Nodes, Layers), 2)

  SupraOutDegree <- sumR(binarizeMatrix_edit(SupraAdjacencyMatrix), 2)
  
  MultiOutDegreeVector <- sumR(reshapeR(SupraOutDegree, Nodes, Layers), 2)

  if (!isDirected) {
    MultiDegreeVector <- (MultiInDegreeVector + MultiOutDegreeVector)/2
  }
  else {
    MultiDegreeVector <- MultiInDegreeVector + MultiOutDegreeVector
  }
  
  return(MultiDegreeVector)
}


binarizeMatrix_edit <- function(SparseMat){
  
  SparseMat@x[SparseMat@x !=  0] <- 1
  
  return(SparseMat)
}

## Eigenvector centrality -------------------------------------------------


GetMultiEigenvectorCentrality_edit <- function(SupraAdjacencyMatrix, Layers, Nodes) {
  
  LeadingEigenvector <- GetLargestEigenv_edit(SupraAdjacencyMatrix, Layers, Nodes, Type = "eigenvector")
  
  CentralityVector <- sumR(reshapeR(LeadingEigenvector, Nodes, Layers), 2)
  
  CentralityVector <- CentralityVector/max(CentralityVector)

  return(CentralityVector)
}

## PageRank centrality ----------------------------------------------------


GetMultiPageRankCentrality_edit <- function(SupraAdjacencyMatrix, Layers, Nodes,
                                            Method = "multilayer") {
  
  LeadingEigenvector <- GetLargestEigenv_edit(SupraAdjacencyMatrix, Layers, Nodes, Type = "pagerank")
  
  CentralityVector <- LeadingEigenvector/sum(LeadingEigenvector)
  
  if (Method == "multilayer") {
    CentralityVector <- sumR(reshapeR(CentralityVector, Nodes, Layers), 2)
  }
  
  CentralityVector <- CentralityVector/max(CentralityVector)
  
  return(CentralityVector)
}


GetLargestEigenv_edit <- function(SupraAdjacencyMatrix, Layers, Nodes, Type = "eigenvector") {
  
  weight <- rep(1, Nodes*Layers)/(Nodes*Layers)
  
  # For pagerank
  if(Type == "pagerank") {
    M_hat <- BuildSupraTransitionMatrixFromSupraAdjacencyMatrix_edit(SupraAdjacencyMatrix)
    M_hat_add <- Matrix(0.85*M_hat, sparse = TRUE)
    extra <- rep(0.15, Nodes*Layers)/(Nodes*Layers)
  # For normal eigenvector
  } else if(Type == "eigenvector"){
    M_hat_add <- SupraAdjacencyMatrix
    extra <- rep(0, Nodes*Layers)
  }
  
  for (i in 1:100) {
    weight <-  (weight %*% M_hat_add)  + (weight %*% extra)[1,1]
  }
  
  return(weight)
}


BuildSupraTransitionMatrixFromSupraAdjacencyMatrix_edit <- function(SupraAdjacencyMatrix) {
  
  degree <- sumR(SupraAdjacencyMatrix, 1)
  D = Diagonal(x = 1/degree)
  A_hat <- D %*% SupraAdjacencyMatrix
  
  return(A_hat)
}


# tau-weighted multi-layer degree centrality ------------------------------


GetMultiTauDegree <- function(SupraAdjacencyMatrix, Layers, Nodes, Taus){
  
  NodesTensor <- SupraAdjacencyToNodesTensor(binarizeMatrix_edit(SupraAdjacencyMatrix), 
                                             Layers, Nodes)
  
  DegreeMatrix <- Matrix(data = 0, nrow = Nodes, ncol = Layers)
  
  for (l in 1:Layers) {
    MultiInDegreeVector <- Matrix::t(sumR(NodesTensor[[l]], 1))
    
    MultiOutDegreeVector <- sumR(NodesTensor[[l]], 2)
    
    DegreeMatrix[,l] <- Taus[l] * ((MultiInDegreeVector + MultiOutDegreeVector)/2)
  }
  
  MultiTauDegree <- Matrix::rowSums(DegreeMatrix)
  
  return(MultiTauDegree)
  
}


# Epidemic modeling -------------------------------------------------------

model_sir_multiplex <- function(adj_mats, tau, gamma, step_size = 7,
                                seed_infs = 1, random_seed = 789, t_max = Inf,
                                verbose = TRUE, out = "time_to inf") {
  # adj_mats = list of separate adj. mats. per layer (= node tensor)
  # tau = vector of probabilities of infection in each layer
  # gamma = recovery time (not layer specific)
  # step_size = length (days) of one iteration
  # seed_infs = no. of infected seed nodes
  # random_seed = seed for sampling seed nodes
  # t_max = msx. no. of time steps
  # verbose = print iterations of while loop?
  # out = return value: sir_mat, time_to_inf, tti_duration, all
  
  layers <- length(adj_mats)
  N <- dim(adj_mats[[1]])[1]
  S <- Matrix(rep(1, N), nrow = N, ncol = 1, sparse = TRUE)
  I <- Matrix(rep(0, N), nrow = N, ncol = 1, sparse = TRUE)
  R <- Matrix(rep(0, N), nrow = N, ncol = 1, sparse = TRUE)
  
  # sample seed infections
  if(!is.null(random_seed)) set.seed(random_seed)
  
  I1 <- sample(1:N, size = seed_infs)
  
  # change status of seed nodes
  S[I1,1] <- 0
  I[I1,1] <- 1
  
  # set recovery time sampled from Weibull distribution per individual
  rec_time <- rweibull(N, scale = gamma, shape = 1)
  
  # storage for new infections per layer
  newI_layers <- vector(mode = "list", length = layers)
  
  # while there are still people infected, i.e., not recovered & t <= t_max...
  t <- 1
  while (sum(I[,t]) > 0 & t <= t_max) {
    
    t <- t + 1
    if(verbose) cat("Step", t, " ")
    
    for (l in 1:layers) {
      #vector of infection status neighbors
      infneigh <- adj_mats[[l]] %*% I[,t-1]
      # prob. of infection along each edge (Reed-Frost model)
      pinf <- as.vector(1 - (1 - tau[l])^infneigh)
      
      # create new infections based on pinf in layer
      newI_layers[[l]] <- rbinom(N, S[,t-1], pinf)
    }
    
    # get union of infections across layers
    # -> infected in one layer = infected in all layers
    newI <- as.integer(apply(do.call("cbind", newI_layers), 1, sum) > 0)
    
    # create new recoveries based on rec_time
    inf_days <- rowSums(I) * step_size
    newR <- as.integer((inf_days - rec_time > 0) & (I[,t-1] == 1))
    
    # update groups for next time step
    nextS <- S[,t-1] - newI
    nextI <- I[,t-1] + newI - newR
    nextR <- R[,t-1] + newR
    
    # collect old and new observations for each group
    S <- cbind(S, nextS)
    I <- cbind(I, nextI)
    R <- cbind(R, nextR)
  }
  
  SIR <- list(S = S, I = I, R = R)
  
  # RETURN based on out argument
  switch (out,
    sir_mat = return(SIR),
    
    epi_curves = {
      cum_S <- unname(apply(SIR$S, 2, sum))
      cum_I <- unname(apply(SIR$I, 2, sum))
      cum_R <- unname(apply(SIR$R, 2, sum))
      return(list("S" = cum_S, "I" = cum_I, "R" = cum_R))},
    
    time_to_inf = {
      # compute times to infection, Inf if not infected
      s_times <- rowSums(S)
      s_times[s_times == ncol(S)] <- Inf
      return(s_times)},
    
    duration = {
      # duration of epidemic
      duration <- ncol(S)
      return(duration)},
    
    tti_duration = {
      s_times <- rowSums(S)
      s_times[s_times == ncol(S)] <- Inf
      duration <- ncol(S)
      return(list("time_to_inf" = s_times, "duration" = duration))},
    
    tti_duration_epi = {
      s_times <- rowSums(S)
      s_times[s_times == ncol(S)] <- Inf
      duration <- ncol(S)
      cum_S <- unname(apply(SIR$S, 2, sum))
      cum_I <- unname(apply(SIR$I, 2, sum))
      cum_R <- unname(apply(SIR$R, 2, sum))
      return(list("time_to_inf" = s_times, "duration" = duration,
                  "S" = cum_S, "I" = cum_I, "R" = cum_R))},
    
    all = {
      s_times <- rowSums(S)
      s_times[s_times == ncol(S)] <- Inf
      duration <- ncol(S)
      return(list("sir_mat" = SIR, "time_to_inf" = s_times, "duration" = duration))}
  )
  
}


# Variable transformations ------------------------------------------------

fisher_z_transform <- function(coef, back = FALSE) {
  # coef = the correlation coefficient to be transformed
  # back = apply back-transformation?
  
  if(!back) {
    coef_z <- 0.5 * log((1 + coef) / (1 - coef))
    return(coef_z)
    
  } else {
    coef_back <- (exp(2 * coef) - 1) / (exp(2 * coef) + 1)
    return(coef_back)
  }
}
  
  

# XGBoost modeling pipeline -----------------------------------------------

xgb_pipeline <- function(seed = 9106, 
                         predictors, outcome, 
                         #sample_size, 
                         eval_metric = "rmse",
                         objective = "reg:squarederror",
                         hyperparams = NULL,
                         train_prop = 0.1) {
  # seed: random seed
  # predictors, outcome: named vectors of variables
  # sample_size: subset of dataset to be used
  # eval_metric: rmse or rsq, e.g., leave out for cox
  # objective: reg:squarederror or survival:cox
  # hyperparams: either NULL to tune within function, or specified hyperparams.
  
  set.seed(seed)
  
  # Data assembly
  xgb_dat <- 
    data.frame(predictors, outcome) %>% 
    drop_na() #%>% 
    #slice_sample(n = sample_size) 
  
  # Data split
  xgb_split <- initial_split(xgb_dat, prop = train_prop)
  xgb_train <- training(xgb_split)
  
  # Hyperparameter specs
  xgb_spec <- boost_tree(
    mode = "regression",
    trees = tune(),
    tree_depth = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    mtry = length(predictors),
    learn_rate = tune()
    ) %>% 
    set_engine("xgboost", 
               objective = objective)
  
  # do tuning if no hyeperparams. are given
  if(is.null(hyperparams)) {
  
    # Tuning grid (space-filling search design)
    xgb_grid <- grid_latin_hypercube(
      trees(range = c(100, 1000)),
      tree_depth(range = c(2,5)),
      min_n(),
      loss_reduction(),
      sample_size = sample_prop(range = c(0.2, 0.9)),
      learn_rate(),
      size = 50
    )
    
    # Workflow
    xgb_wf <- workflow() %>% 
      add_formula(tti ~.) %>% 
      add_model(xgb_spec)
    
    # Create CV folds
    xgb_folds <- vfold_cv(xgb_train)
    
    # Do the tuning
    all_cores <- parallel::detectCores()
    
    cl <- makePSOCKcluster(all_cores)
    registerDoParallel(cl)
    
    xgb_res <- tune_grid(
      xgb_wf,
      resamples = xgb_folds,
      grid = xgb_grid,
      control = control_grid(save_pred = FALSE)
    )
  
    stopCluster(cl)
    
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
    
    # Select model with best hyperparameter combination
    best_params <- select_best(xgb_res, metric = eval_metric)

  } else {
    
    xgb_wf <- workflow() %>% 
      add_formula(tti ~.) %>% 
      add_model(xgb_spec)
    
    best_params <- hyperparams
  }
  
  # Add hyperparams. to workflow
  xgb_final <- finalize_workflow(xgb_wf, best_params)
  
  # Fit best model on training data and evaluate on test data
  xgb_res_final <- last_fit(xgb_final, xgb_split)
  
  # Get performance metrics
  metrics_final <- collect_metrics(xgb_res_final)
  
  # Get variable importance measures
  # (no work when data too small)
  vi <- tryCatch(xgb_res_final %>% 
                   extract_workflow() %>%
                   extract_fit_engine() %>%
                   xgb.importance(model = .),
                  error = function(e) {
                    paste(e$message)
                    FALSE},
                  warning = function(w) {
                    paste(w$message)
                    FALSE
                    })
  
  return(list("model_params" = best_params,
              "metrics" = metrics_final, 
              "vimportance" = vi,
              "last_fit" = extract_workflow(xgb_res_final)))
}


  
  
  
  
  
  
  





























  