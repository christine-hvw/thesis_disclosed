# _____________________________________________________________________________
# title: "Predicting Covid-19 infections using multi-layer centrality measures"
# subtitle: "Auxiliary functions"
# author: "Christine Hedde-von Westernhagen"
# date: "15-05-2023"
# _____________________________________________________________________________


# Testing -----------------------------------------------------------------
# # >>> UNCOMMENT THIS PART FOR TESTING OF FUNCTIONS <<<
# Packages
# library(muxViz)
# library(Matrix)
# library(igraph)
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

# All code based on the muxViz package (v3.1) as stated below has been acquired via 
# https://github.com/manlius/muxViz, which is licensed and copyrighted under 
# GNU General Public License v3.0.

## Degree centrality ------------------------------------------------------

# Based on muxViz::GetMultiDegree()
# - combines  GetMultiInDegree() and  GetMultiOutDegree(), uses adaption of binarizeMatrix()

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


## Binarize matrix --------------------------------------------------------

# Based on muxViz::binarizeMatrix()
# - adaption necessary because it throws error with large matrices

binarizeMatrix_edit <- function(SparseMat){
  
  SparseMat@x[SparseMat@x !=  0] <- 1
  
  return(SparseMat)
}

## Eigenvector centrality -------------------------------------------------

# Based on muxViz::GetMultiEigenvectorCentrality()
# - uses power method from GetLargestEigenv_edit() to find eigenvector centrality

GetMultiEigenvectorCentrality_edit <- function(SupraAdjacencyMatrix, Layers, Nodes) {
  
  LeadingEigenvector <- GetLargestEigenv_edit(SupraAdjacencyMatrix, Layers, Nodes, 
                                              Type = "eigenvector")
  
  CentralityVector <- sumR(reshapeR(LeadingEigenvector, Nodes, Layers), 2)
  
  CentralityVector <- CentralityVector/max(CentralityVector)

  return(CentralityVector)
}

## PageRank centrality ----------------------------------------------------

# Based on muxViz::GetMultiEigenvectorCentrality()
# - uses power method from GetLargestEigenv_edit() to find PageRank centrality

GetMultiPageRankCentrality_edit <- function(SupraAdjacencyMatrix, Layers, Nodes,
                                            Method = "multilayer") {
  
  LeadingEigenvector <- GetLargestEigenv_edit(SupraAdjacencyMatrix, Layers, Nodes, 
                                              Type = "pagerank")
  
  CentralityVector <- LeadingEigenvector/sum(LeadingEigenvector)
  
  if (Method == "multilayer") {
    CentralityVector <- sumR(reshapeR(CentralityVector, Nodes, Layers), 2)
  }
  
  CentralityVector <- CentralityVector/max(CentralityVector)
  
  return(CentralityVector)
}


## Get largest eigenvector ------------------------------------------------

# Replaces muxViz::GetLargestEigenv() in Eigenv. and PageRank centrality functions
# - power method is used which facilitates finding eigenvectors in large matrices

GetLargestEigenv_edit <- function(SupraAdjacencyMatrix, Layers, Nodes, Type = "eigenvector") {
  
  # Initial guess of (normalized) eigenvector
  weight <- rep(1, Nodes * Layers)/(Nodes * Layers)
  
  # For PageRank centrality
  if(Type == "pagerank") {
    
    # Get transition matrix of supra adjacency matrix
    M_hat <- BuildSupraTransitionMatrixFromSupraAdjacencyMatrix_edit(SupraAdjacencyMatrix)
    
    # Add PageRank factor
    M_hat_add <- Matrix(0.85 * M_hat, sparse = TRUE)
    
    extra <- rep(0.15, Nodes * Layers)/(Nodes * Layers)
    
  # For normal eigenvector don't add anything
  } else if(Type == "eigenvector"){
    
    M_hat_add <- SupraAdjacencyMatrix
    
    extra <- rep(0, Nodes*Layers)
  }
  
  # Repeatedly multiply eigenvector guess with adjacency matrix
  # (using commutative property)
  for (i in 1:100) {
    
    weight <-  (weight %*% M_hat_add)  + (weight %*% extra)[1,1]
  }
  
  return(weight)
}


## Supra-transition matrix ------------------------------------------------

# Reduced version of muxViz::BuildSupraTransitionMatrixFromSupraAdjacencyMatrix()
# - used in combination with adapted functions for Eigenv. and PageRank centrality

BuildSupraTransitionMatrixFromSupraAdjacencyMatrix_edit <- function(SupraAdjacencyMatrix) {
  
  degree <- sumR(SupraAdjacencyMatrix, 1)
  D = Diagonal(x = 1/degree)
  A_hat <- D %*% SupraAdjacencyMatrix
  
  return(A_hat)
}


## Tau-weighted multi-layer Degree centrality ------------------------------

# Version of multi-layer Degree centrality with layer specific weighting ("Taus")

GetMultiTauDegree <- function(SupraAdjacencyMatrix, Layers, Nodes, Taus){
  
  NodesTensor <- SupraAdjacencyToNodesTensor(binarizeMatrix_edit(SupraAdjacencyMatrix), 
                                             Layers, Nodes)
  
  # Create storage
  DegreeMatrix <- Matrix(data = 0, nrow = Nodes, ncol = Layers)
  
  for (l in 1:Layers) {
    # Compute in- and out-degree
    MultiInDegreeVector <- Matrix::t(sumR(NodesTensor[[l]], 1))
    
    MultiOutDegreeVector <- sumR(NodesTensor[[l]], 2)
    
    # Weight Degree by layer
    DegreeMatrix[,l] <- Taus[l] * ((MultiInDegreeVector + MultiOutDegreeVector)/2)
  }
  
  # Get sum of weighted Degrees
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
  
  # Set up parameters and storage for infection status groups
  layers <- length(adj_mats)
  N <- dim(adj_mats[[1]])[1]
  S <- Matrix(rep(1, N), nrow = N, ncol = 1, sparse = TRUE)
  I <- Matrix(rep(0, N), nrow = N, ncol = 1, sparse = TRUE)
  R <- Matrix(rep(0, N), nrow = N, ncol = 1, sparse = TRUE)
  
  # Sample seed infections
  if(!is.null(random_seed)) set.seed(random_seed)
  
  I1 <- sample(1:N, size = seed_infs)
  
  # Change status of seed nodes
  S[I1,1] <- 0
  I[I1,1] <- 1
  
  # Set recovery time sampled from Weibull distribution per individual
  rec_time <- rweibull(N, scale = gamma, shape = 1)
  
  # Storage for new infections per layer
  newI_layers <- vector(mode = "list", length = layers)
  
  # While there are still people infected, i.e., not recovered & t <= t_max...
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
    # (3 matrices of N x t)
    S <- cbind(S, nextS)
    I <- cbind(I, nextI)
    R <- cbind(R, nextR)
  }
  
  # Collect all status group matrices
  SIR <- list(S = S, I = I, R = R)
  
  # RETURN based on "out" argument
  switch (out,
    
    # complete matrices of infection status over time
    sir_mat = return(SIR),
    
    # times to infection, Inf if not infected
    time_to_inf = {
      s_times <- rowSums(S)
      s_times[s_times == ncol(S)] <- Inf
      return(s_times)},
    
    # duration of epidemic
    duration = {
      duration <- ncol(S)
      return(duration)},
    
    # both tti and duration
    tti_duration = {
      s_times <- rowSums(S)
      s_times[s_times == ncol(S)] <- Inf
      duration <- ncol(S)
      return(list("time_to_inf" = s_times, "duration" = duration))},
    
    # all of the above
    all = {
      s_times <- rowSums(S)
      s_times[s_times == ncol(S)] <- Inf
      duration <- ncol(S)
      return(list("sir_mat" = SIR, "time_to_inf" = s_times, "duration" = duration))}
  )
  
}


# Variable transformations ------------------------------------------------

# Perform a Fisher z-transformation to normality on a variable (or back-transformation)

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
  
  



  
  
  
  
  
  
  





























  