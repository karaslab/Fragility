require(glmnet)
require(doMC)

load_fragility_patient <- function(subject_code, block, elec, unit, halve = FALSE) {
  v2 <- readRDS(paste0('/Volumes/bigbrain/oliver-r-projects/', subject_code, ' R Data/', subject_code, '_car_voltage'))
  preload_info <- readRDS(paste0('/Volumes/bigbrain/oliver-r-projects/', subject_code, ' R Data/', subject_code, '_car_info'))
  
  cond = preload_info$condition
  
  elec <- as.character(elec)
  trial_num <- match(block,cond)
  # trial_num <- match(paste0('seizure (',block,')'), cond)
  # N <- length(elec) # number of electrodes
  
  v1 <- v2[,,elec] # assign only specified electrodes
  
  if (unit != 'uV') {
    if (unit == 'nV') {
      v1 <- v1/1000
    } else if (unit == 'mV') {
      v1 <- v1*1000
    } else {
      stop('Accepted units are uV, nV, or mV')
    }
  }
  
  
  if (halve) {
    v1 <- v1[,seq(1,ncol(v2),2),] # halves the frequency
  }
  
  pt_info <- list(
    v = v1,
    trial = trial_num
  )
}

generate_adj_array <- function(t_window, t_step, v, trial_num, nlambda, ncores) {
  S <- dim(v)[2] # S is total number of timepoints
  N <- dim(v)[3] # N is number of electrodes
  
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
    svec <- generate_state_vectors(v,trial_num,t_window,t_start)
    A[,,k] <- find_adj_matrix(svec, N, t_window, nlambda = nlambda, ncores = ncores)
    
    # MSE
    estimate <- A[,,k] %*% svec$x
    mse[k] <- mean((estimate - svec$x_n)^2)
    
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  
  adj_info <- list(
    A = A,
    mse = mse
  )
}

generate_fragility_matrix <- function(A, elec, lim = 1i) {
  N <- dim(A)[1]
  J <- dim(A)[3]
  f_vals <- matrix(nrow = N, ncol = J)
  for (k in 1:J) {
    start_time <- Sys.time()
    print(paste0('current timewindow: ', k, ' out of ', J))
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

find_adj_matrix <- function(state_vectors, N, t_window, nlambda, ncores) {
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
  cv.ridge <- glmnet::cv.glmnet(H, b, alpha = 0, nfolds = 3, parallel = TRUE, nlambda = nlambda)
  lambdas <- rev(cv.ridge$lambda)
  
  test_lambda <- function(l, H, b, lambdas, ncores) {
    ridge <- glmnet::glmnet(H, b, alpha = 0, lambda = l)
    N = sqrt(dim(H)[2])
    adj_matrix <-  matrix(ridge$beta, nrow = N, ncol = N, byrow = TRUE)
    eigv <- abs(eigen(adj_matrix, only.values = TRUE)$values)
    stable <- max(eigv) < 1
    list(
      adj = adj_matrix,
      abs_eigv = eigv,
      stable = stable
    )
  }
  
  l <- 1
  stable_i <- vector(length = 0)
  
  while (length(stable_i) == 0) {
    
    if ((l+ncores-1) <= length(lambdas)) {
      results <- rave::lapply_async3(lambdas[l:(l+ncores-1)], test_lambda, H = H, b = b, .ncores = ncores)
      stable_i <- which(unname(unlist(lapply(results, function (x) x['stable']))))
    } else {
      results <- rave::lapply_async3(lambdas[l:length(lambdas)], test_lambda, H = H, b = b, .ncores = ncores)
      stable_i <- which(unname(unlist(lapply(results, function (x) x['stable']))))
      break
    }
    
    l <- l + ncores
  }
  
  if (length(stable_i) == 0) {
    stop('no lambdas make result in stable adjacency matrix')
  }
  
  adj_matrix <- results[[stable_i[1]]]$adj
  
  return(adj_matrix)
}

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

requested_electrodes <- c(1:24,26:36,42:43,46:54,56:70,72:95)

pt_info <- load_fragility_patient(
  subject_code = 'PT01',
  block = 'seizure (1)',
  elec = requested_electrodes,
  unit = 'nV'
)

adj_info <- generate_adj_array(
  t_window = 250, 
  t_step = 250, 
  v = pt_info$v, 
  trial_num = pt_info$trial,
  nlambda = 16,
  ncores = 8
)

f_info <- generate_fragility_matrix(
  A = adj_info$A,
  elec = requested_electrodes
)

f_plot_params <- list(
  mat = f_info$norm,
  x = 1:J,
  y = 1:N,
  zlim = c(0,1)
)

ridge <- glmnet::glmnet(H, b, alpha = 0, lambda = 50)
adj_matrix <-  matrix(ridge$beta, nrow = N, ncol = N, byrow = TRUE)
max(abs(eigen(adj_matrix, only.values = TRUE)$values))

image(z = t(f_info$norm), x=1:160,
      y = 1:N, # as.numeric(elec),
      zlim=c(0,1),
      col = colorRampPalette(c('blue', 'green', 'red'))(101))
