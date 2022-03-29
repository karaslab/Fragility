require(glmnet)
require(doMC)

# functions ----

# generate state vectors
# x(t) represents the voltages of each electrode at time t in a N by T-1 matrix
# where N is # of electrodes and T is timepoints in the specified time window
# x(t+1) represents the same thing except shifted over by one timepoint
# t_start specifies starting timepoint for this time window
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
  #state_vectors <- svec # for testing purposes
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
  cv.ridge <- cv.glmnet(H, b, alpha = 0, nfolds = 3, parallel = TRUE)
  l <- length(cv.ridge$lambda) + 1
  
  eigv <- 1:63
  
  # solve least squares for eigenvalues of adj_matrix less than 1
  while (max(eigv) >= 1) {
    
    if (l == 1) {
      stop('no lambdas make abs(eigenvalues) less than 1')
    }
    
    l <- l - 1
    lambda <- cv.ridge$lambda[l]
    
    ridge <- glmnet(H, b, alpha = 0, lambda = lambda)
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

# lower output value means higher fragility

# load pt info ----
# # YDR
# v <- readRDS('/Volumes/bigbrain/oliver-r-projects/YDR R Data/YDR_car_voltage_all')
# preload_info <- readRDS('/Volumes/bigbrain/oliver-r-projects/YDR R Data/YDR_car_info_all')
# 
# cond <- preload_info$condition
# block <- 'V'
# elec <- as.character(c(61:87,93:152))

# YDS grid
# v2 <- readRDS('~/Documents/R/YDS_car_voltage_grid')
# preload_info <- readRDS('~/Documents/R/YDS_car_info_grid')
# 
# cond <- preload_info$condition
# block <- 'B'
# elec <- as.character(132:194)
# electrodes 174:180 result in a stable adjacency matrix

# YDS seeg
v2 <- readRDS('~/Documents/R/YDS_car_voltage_seeg')
preload_info <- readRDS('~/Documents/R/YDS_car_info_seeg')

cond <- preload_info$condition
block <- 'B'
elec <- as.character(c(1:19,21:35,42:131,196:207))

# PT01
# v2 <- readRDS('~/Documents/R/PT01_car_voltage')
# v2 <- v2/1000 # convert from nanovolts to microvolts
# preload_info <- readRDS('~/Documents/R/PT01_car_info')
# 
# cond <- preload_info$condition
# block <- '1'
# elec <- as.character(c(1:24,26:36,42:43,46:54,56:70,72:95))

# KAA
# v2 <- readRDS('~/Documents/R/KAA_car_voltage')
# v2 <- v2 * 1000 # convert from millivolts to microvolts
# preload_info <- readRDS('~/Documents/R/KAA_car_info')

# cond <- preload_info$condition
# block <- '4'
# elec <- as.character(c(1:43,45:80,82:116,129:164,166:244))
# elec <- as.character(c(1:43,45,46,61:72,89:116,145:158,196:209))

trial_num <- match(paste0('seizure (',block,')'), cond)
N <- length(elec) # number of electrodes

v1 <- v2[,,elec] # assign only specified electrodes
v1 <- v1[,seq(1,ncol(v2),2),] # if 2000Hz, halve to 1000Hz

# code ----
# length of a single time window in whichever units one timepoint is
t_window <- 250
t_step <- 125

# find adjacency matrices for entire timecourse, one for each time window
S <- ncol(v1) # S is total number of timepoints
if(S %% t_step != 0) {
  # truncate S to greatest number evenly divisible by timestep
  S <- trunc(S/t_step) * t_step
}
J <- S/t_step - (t_window/t_step) + 1 # J is number of time windows
A <- array(dim = c(N,N,J))
mse <- vector(mode = "numeric", length = J)

ncores = 4
registerDoMC(cores = ncores)

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
mean(mse)

abs(eigen(A[,,1])$values)

saveRDS(A, file = '/Volumes/bigbrain/oliver-r-projects/KAA_adj_laptop')
A_s <- readRDS('/Volumes/bigbrain/oliver-r-projects/KAA R Data/KAA_adj_server')

# reconstruct model using A to see if stable?
v_trace <- t(v1[trial_num,,])
v_recon <- v_trace
for (k in 1:J){
  t_start <- 1+(k-1)*t_step
  v_recon[,t_start] <- v_trace[,t_start]
  for (i in (t_start+1):(t_start+t_window-1)) {
    v_recon[,i+1] <- A[,,k] %*% v_recon[,i]
  }
}
for (i in 1:(S-1)) {
  v_recon[,i+1] <- A[,,1] %*% v_recon[,i]
}

max(abs(eigen(A[,,1], only.values = TRUE)$values))

timepoints <- 1:S
plot(x = timepoints, y = v_trace[1,timepoints], type = 'l')
plot(x = timepoints, y = v_recon[1,timepoints], type = 'l')

# find fragility of A
f_vals <- matrix(nrow = N, ncol = J)
lim <- 1i
for (k in 1:J) {
  start_time <- Sys.time()
  print(paste0('current index: ', k, ' out of ', J))
  for (i in 1:N) {
    f_vals[i,k] <- find_fragility(i,A[,,k],N,lim)
  }
  end_time <- Sys.time()
  print(end_time - start_time)
}
rownames(f_vals) <- attr(v1, 'dimnames')$Electrode

saveRDS(f_vals, file = '/Volumes/bigbrain/oliver-r-projects/KAA_f_vals')

f_normalized <- f_vals

# scale fragility values from 0 to 1 with 1 being most fragile
for (j in 1:J) {
  max_f <- max(f_vals[,j])
  f_normalized[,j] <- sapply(f_vals[,j], function(x) (max_f - x) / max_f)
}

image(z = t(f_normalized), x=1:160,
      y = 1:N, # as.numeric(elec),
      zlim=c(0,1),
      col = colorRampPalette(c('blue', 'green', 'red'))(101))

avg_f <- rowMeans(f_normalized)

# export avg_f
# df <- data.frame(avg_f)
# write.csv(df,'~/Documents/R/f_avg.csv')

# test space----

f_test <- f_normalized[as.character(c(34:35,63:70,89:92)),]
f_flip <- f_test[nrow(f_test):1,]
image(z = t(f_flip), x=1:160,
      y = 1:14,
      zlim=c(0,1),
      col = colorRampPalette(c('blue', 'green', 'red'))(101))
