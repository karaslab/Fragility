#' load_fragility_patient
#' this function loads a previously processed patient
#' @param subject_code Subject code
#' @param block The epoch condition/trial that you're calculating for
#' @param elec Vector of electrodes you want included
#' @param halve Logical determining whether you want to halve the Hz value, e.g. from 2000Hz to 1000Hz (saves some computation time)
#'
#' @return list containing the voltage trace matrix v, condition,
#' @export
#'
#' @examples
load_fragility_patient <- function(subject_code, block, elec, halve = FALSE) {
  v2 <- readRDS(paste0('/Volumes/bigbrain/oliver-r-projects/', subject_code, ' R Data/', subject_code, '_car_voltage'))
  preload_info <- readRDS(paste0('/Volumes/bigbrain/oliver-r-projects/', subject_code, ' R Data/', subject_code, '_car_info'))
  
  elec <- as.character(elec)
  trial_num <- match(block,cond)
  # trial_num <- match(paste0('seizure (',block,')'), cond)
  N <- length(elec) # number of electrodes
  
  v1 <- v2[,,elec] # assign only specified electrodes
  
  if (halve) {
    v1 <- v1[,seq(1,ncol(v2),2),] # halves the frequency
  }
  
  pt_info <- list(
    v = v1,
    cond = preload_info$condition,
    trial = trial_num,
    N = N
  )
}

generate_adj_array <- function(t_window, t_step, v, ncores) {
  S <- ncol(v) # S is total number of timepoints
  if(S %% t_step != 0) {
    # truncate S to greatest number evenly divisible by timestep
    S <- trunc(S/t_step) * t_step
  }
  J <- S/t_step - (t_window/t_step) + 1 # J is number of time windows
  A <- array(dim = c(N,N,J))
  mse <- vector(mode = "numeric", length = J)
  
  doMC::registerDoMC(cores = ncores)
  
  for (k in 1:J) {
    start_time <- Sys.time()
    print(paste0('current timewindow: ', k, ' out of ', J))
    t_start <- 1+(k-1)*t_step
    svec <- generate_state_vectors(v1,trial_num,t_window,t_start)
    A[,,k] <- find_adj_matrix(svec, N, t_window)
    
    # MSE
    estimate <- A[,,k] %*% svec$x
    mse[k] <- mean((estimate - svec$x_n)^2)
    
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  
  adj_info <- list(
    A = A,
    mse = mse,
    J = J
  )
}

generate_fragility_matrix <- function(N, J, A, elec, lim = 1i) {
  f_vals <- matrix(nrow = N, ncol = J)
  for (k in 1:J) {
    start_time <- Sys.time()
    print(paste0('current index: ', k, ' out of ', J))
    for (i in 1:N) {
      f_vals[i,k] <- find_fragility(i,A[,,k],N,lim)
    }
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  rownames(f_vals) <- elec
  colnames(f_vals) <- 1:J
  
  f_norm <- f_vals
  
  # scale fragility values from 0 to 1 with 1 being most fragile
  for (j in 1:J) {
    max_f <- max(f_vals[,j])
    f_norm[,j] <- sapply(f_vals[,j], function(x) (max_f - x) / max_f)
  }
  
  avg_f <- rowMeans(f_norm)
  
  f_info <- list(
    vals = f_vals,
    norm = f_norm,
    avg = avg_f
  )
}

generate_state_vectors <- function(v,trial,t_window,t_start) {
  data <- v[trial,,]
  
  state_vectors <- list(
    # x(t)
    x = t(data)[,t_start:(t_start+t_window-2)],
    
    # x(t+1)
    x_n = t(data)[,(t_start+1):(t_start+t_window-1)]
  )
  return(state_vectors)
}

# finds the adjacency matrix for given time window
find_adj_matrix <- function(state_vectors, N, t_window) {
  # vectorize x(t+1)
  # state_vectors <- svec # for testing purposes
  b <- c(state_vectors$x_n)
  
  # initialize big H matrix for system of linear equations
  H <- matrix(0, nrow = N*(t_window-1), ncol = N^2)
  
  # populate H matrix
  r <- 1
  for (ii in 1:(t_window-1)) {
    c <- 1
    for (jj in 1:N) {
      H[r,c:(c+N-1)] <- state_vectors$x[,ii]
      c <- c + N
      r <- r + 1
    }
  }
  
  # solve system using glmnet package least squares, with L2-norm regularization
  # aka ridge filtering
  
  # find optimal lambda
  cv.ridge <- glmnet::cv.glmnet(H, b, alpha = 0, nfolds = 3, parallel = TRUE)
  l <- length(cv.ridge$lambda) + 1
  
  eigv <- 1:63
  
  # solve least squares for eigenvalues of adj_matrix less than 1
  while (max(eigv) >= 1) {
    
    if (l == 1) {
      stop('no lambdas make abs(eigenvalues) less than 1')
    }
    
    l <- l - 1
    lambda <- cv.ridge$lambda[l]
    
    ridge <- glmnet::glmnet(H, b, alpha = 0, lambda = lambda)
    X <- ridge$beta
    
    adj_matrix <- matrix(X, nrow = N, ncol = N, byrow = TRUE)
    eigv <- abs(eigen(adj_matrix, only.values = TRUE)$values)
  }
  # X <- solve(t(H) %*% H) %*% t(H) %*% b
  # adj_matrix <- matrix(X, nrow = N, ncol = N, byrow = TRUE)
  
  return(adj_matrix)
}

#' find_fragility
#' Finds the fragility value for a single electrode (node) in a single time window
#'
#' @param node The node that you're calculating the fragility for
#' @param A The adjacency matrix that you're calculating the fragility for
#' @param N The number of electrodes
#'
#' @return
#' @export
#'
#' @examples

find_fragility <- function(node, A_k, N, limit) {
  
  e_k <- vector(mode = 'numeric', length = N)
  e_k[node] <- 1
  # limit <- 0.707107+0.707107i
  # limit <- 1+1i
  
  argument <- t(e_k) %*% (solve(A_k - limit*diag(N))) # column perturbation
  # argument <- t(e_k) %*% t(solve(A_k - num*diag(N))) # row perturbation
  
  B <- rbind(Im(argument),Re(argument))
  
  perturb_mat <- (t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) %*% t(e_k) # column
  # perturb_mat <- e_k %*% t(t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) # row
  
  norm(perturb_mat, type = '2')
}