#' process_fragility_patient
#'
#' This function processes voltage data in preparation for adjacency matrix and
#' fragility calculations. It converts the units of the data to uV if possible
#' and halves the frequency if necessary.
#'
#' @param v 3D voltage array from RAVE's module_tools$get_voltage() and
#'   get_data() functions. First dimension is trial, second dimension is time,
#'   and third dimension is electrode.
#' @param unit Character: 'uV', 'nV', or 'mV' specifying the units of the
#'   voltage data.
#' @param halve Logical specifying whether the data should be halved or not.
#'   This is for data that has a high sampling rate - for example, 2000 Hz data
#'   would be processed much more efficiently if halved to 1000Hz, while still
#'   retaining most of the important information.
#'
#' @return List containing the processed voltage trace matrix v as well as
#'   metadata about the PROCESSED data (trial numbers, unit, sampling rate) -
#'   NOT the original data!
#' @export
#'
#' @examples
#' voltage <- module_tools$get_voltage()
#' v <- voltage$get_data()
#' pt_info <- process_fragility_patient(v = v, unit = 'mV', halve = FALSE)
process_fragility_patient <- function(v, unit, srate, halve = FALSE) {
  print('Pre-processing fragility patient')
  
  newunit <- 'uV'
  
  if (unit != 'uV') {
    if (unit == 'nV') {
      v <- v/1000
    } else if (unit == 'mV') {
      v <- v*1000
    } else {
      warning('Accepted units are uV, nV, or mV')
      newunit <- unit
    }
  }
  
  if (halve) {
    v <- v[,seq(1,ncol(v),2),] # halves the frequency
    srate <- srate/2
  }
  
  pt_info <- list(
    v = v,
    trials = as.numeric(attr(v, "dimnames")[[1]]),
    unit = newunit,
    srate = srate
  )
}

#' generate_adj_array
#'
#' This function generates the NxNxJ adjacency array, where N is number of
#' electrodes and J is number of timewindows. Each NxN matrix within the array
#' is a linear least-squares approximation of how the voltage data evolves
#' within a specific timewindow J. In its entirety, this array can be used as a
#' linear time-varying model of the EEG data, and is used for creating the
#' fragility map.
#'
#' @param t_window Integer specifying size of one time window. Units are in
#'   whatever increment the timepoints are - for example, for 1000Hz data, units
#'   would be in ms.
#' @param t_step Integer specifying size of one time step (when t_step =
#'   t_window, adjacent time windows have no overlap)
#' @param v 3D voltage array from RAVE's module_tools$get_voltage() and
#'   get_data() functions after processing from process_fragility_patient. First
#'   dimension is trial, second dimension is time, and third dimension is
#'   electrode.
#' @param trial_num Integer specifying which trial to create ajdacency array
#'   for. The adjacency array can only be produced with one trial at a time.
#' @param nlambda Integer specifying how many lambdas are to be calculated
#'   during cross-validation for L2-norm regularization (Ridge regression). More
#'   lambdas results in better fitting of data, at the expense of processing
#'   time. Default is 16 lambdas.
#' @param ncores Integer specifying how many cores to utilize during parallel
#'   processing.
#'
#' @return List containing the adjacency array A as well as a measure of how
#'   accurate the fitted data is, mse
#' @export
#'
#' @examples
#' adj_info <- generate_adj_array(
#'     t_window = 250,
#'     t_step = 125,
#'     v = pt_info$v,
#'     trial_num = 1,
#'     nlambda = 16,
#'     ncores = 8
#' )
generate_adj_array <- function(t_window, t_step, v, trial_num, nlambda, ncores) {
  print('Generating adjacency array')
  
  S <- dim(v)[2] # S is total number of timepoints
  N <- dim(v)[3] # N is number of electrodes
  
  if(S %% t_step != 0) {
    # truncate S to greatest number evenly divisible by timestep
    S <- trunc(S/t_step) * t_step
  }
  J <- S/t_step - (t_window/t_step) + 1 # J is number of time windows
  A <- array(dim = c(N,N,J))
  mse <- vector(mode = "numeric", length = J)
  
  # doMC::registerDoMC(cores = ncores)
  
  adjprogress = rave::progress(title = 'Generating Adjacency Array', max = J)
  shiny::showNotification('Calculating estimated time remaining...', id = 'first_est', duration = NULL)
  
  for (k in seq(1,J,ncores)) {
    if (k+ncores-1 <= J) {
      ks <- k:(k+ncores-1)
      
      start_time <- Sys.time()
      if (ncores == 1) {
        print(paste0('Current timewindow: ', k, ' out of ', J))
        adjprogress$inc(paste0('Current timewindow: ', k, ' out of ', J))
      } else {
        print(paste0('Current timewindows: ', k, '-', k+ncores-1, ' out of ', J))
        for (i in ks) {
          adjprogress$inc(paste0('Current timewindows: ', k, '-', k+ncores-1, ' out of ', J))
        }
      }
      
      t_start <- 1+(ks-1)*t_step
      # no significant diff between parallel and sequential for this calculation
      # svec <- lapply(t_start, generate_state_vectors, v = v, trial = trial_num, t_window = t_window)
      svec <- rave::lapply_async3(t_start, generate_state_vectors, v = v, trial = trial_num, t_window = t_window, .ncores = ncores)
      
      A_list <- rave::lapply_async3(svec, find_adj_matrix, N = N, t_window = t_window, nlambda = nlambda, .ncores = ncores)
      
      A[,,ks] <- array(unlist(A_list), dim = c(N,N,length(ks)))
      
      # MSE
      for (i in 1:length(ks)) {
        estimate <- A[,,ks[i]] %*% svec[[i]]$x
        mse[ks[i]] <- mean((estimate - svec[[i]]$x_n)^2)
      }
      
      end_time <- Sys.time()
      print(end_time - start_time)
      
      if (k == 1) {
        shiny::removeNotification(id = 'first_est')
        t_avg <- 0
      }
      
      t_avg <- (t_avg*(((k-1)/ncores)) + as.numeric(difftime(end_time, start_time, units='mins')))/(((k-1)/ncores)+1)
      shiny::showNotification(paste0('Estimated time remaining: ', ((t_avg/ncores)*(J-k))%/%60, ' hours, ', round(((t_avg/ncores)*(J-k))%%60, digits = 1), ' minutes'), id = 'est_time', duration = NULL)
      
    } else {
      ks <- k:J
      
      start_time <- Sys.time()
      print(paste0('Current timewindows: ', k, '-', J, ' out of ', J))
      for (i in ks) {
        adjprogress$inc(paste0('Current timewindows: ', k, '-', k+ncores-1, ' out of ', J))
      }
      
      t_start <- 1+(ks-1)*t_step
      # no significant diff between parallel and sequential for this calculation
      # svec <- lapply(t_start, generate_state_vectors, v = v, trial = trial_num, t_window = t_window)
      svec <- rave::lapply_async3(t_start, generate_state_vectors, v = v, trial = trial_num, t_window = t_window, .ncores = ncores)
      
      A_list <- rave::lapply_async3(svec, find_adj_matrix, N = N, t_window = t_window, nlambda = nlambda, ncores = ncores, .ncores = ncores)
      
      A[,,ks] <- array(unlist(A_list), dim = c(N,N,length(ks)))
      
      # MSE
      for (i in 1:length(ks)) {
        estimate <- A[,,ks[i]] %*% svec[[i]]$x
        mse[ks[i]] <- mean((estimate - svec[[i]]$x_n)^2)
      }
      
      end_time <- Sys.time()
      print(end_time - start_time)
    }
  }
  
  shiny::removeNotification(id = 'est_time')
  
  adjprogress$close()
  
  adj_info <- list(
    A = A,
    mse = mse
  )
}

#' generate_fragility_matrix
#'
#' This function generates the NxJ fragility matrix, where N is the number of
#' electrodes and J is the number of timewindows.
#'
#' @param A NxNxJ adjacency array generated by generate_adj_array.
#' @param elec Character or integer vector containing all the included
#'   electrodes.
#' @param lim Optional parameter setting the eigenvalue used for fragility
#'   calculations. Must be an imaginary or complex number with a absolute value
#'   greater than or equal to 1. Default lim is 1i.
#'
#' @return List containing the raw fragility data (lower number = more fragile),
#'   normalized fragility data from 0 to 1 (1 is most fragile), and average
#'   fragility data across all time windows.
#' @export
#'
#' @examples
#' f_info <- generate_fragility_matrix(
#'     A = adj_info$A, 
#'     elec = c(1:24,26:36,42:43,46:54,56:70,72:95)
#' )
generate_fragility_matrix <- function(A, elec, lim = 1i) {
  print('Generating fragility matrix')
  
  N <- dim(A)[1]
  J <- dim(A)[3]
  f_vals <- matrix(nrow = N, ncol = J)
  fprogress = rave::progress(title = 'Generating Fragility Matrix', max = J)
  shiny::showNotification('Calculating estimated time remaining...', id = 'first_est', duration = NULL)
  for (k in 1:J) {
    start_time <- Sys.time()
    print(paste0('Current timewindow: ', k, ' out of ', J))
    fprogress$inc(paste0('Current timewindow: ', k, ' out of ', J))
    for (i in 1:N) {
      f_vals[i,k] <- find_fragility(i,A[,,k],N,lim)
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    
    if (k == 1) {
      shiny::removeNotification(id = 'first_est')
      t_avg <- 0
    }
    
    t_avg <- (t_avg*(k-1) + as.numeric(difftime(end_time, start_time, units='sec')))/k
    shiny::showNotification(paste0('Estimated time remaining: ', (t_avg*(J-k))%/%60, ' minutes'), id = 'est_time', duration = NULL)
    
  }
  shiny::removeNotification(id = 'est_time')
  fprogress$close()
  rownames(f_vals) <- elec
  colnames(f_vals) <- 1:J
  
  f_norm <- f_vals
  
  # # scale fragility values from 0 to 1 with 1 being most fragile
  # for (j in 1:J) {
  #   max_f <- max(f_vals[,j])
  #   f_norm[,j] <- sapply(f_vals[,j], function(x) (max_f - x) / max_f)
  # }
  
  # scale fragility values from -1 to 1 with 1 being most fragile
  for (j in 1:J) {
    max_f <- max(f_vals[,j])
    f_norm[,j] <- sapply(f_vals[,j], function(x) 2*(max_f - x)/max_f - 1)
  }
  
  avg_f <- rowMeans(f_norm)
  
  f_info <- list(
    vals = f_vals,
    norm = f_norm,
    avg = avg_f
  )
}

#' generate_state_vectors
#'
#' Generates state vectors x(t) and x(t+1) for use in finding the adjacency
#' matrix A in the equation x(t+1) = Ax(t). For use within the
#' generate_adj_array function.
#'
#' @param v 3D voltage matrix
#' @param trial trial number
#' @param t_window timewindow size
#' @param t_start start time of current timewindow
#'
#' @return
#'
#' @examples
#' state_vectors <- generate_state_vectors(v,trial = 1,t_window = 250,t_start = 1)
generate_state_vectors <- function(t_start,v,trial,t_window) {
  data <- v[trial,,]
  
  state_vectors <- list(
    # x(t)
    x = t(data)[,t_start:(t_start+t_window-2)],
    
    # x(t+1)
    x_n = t(data)[,(t_start+1):(t_start+t_window-1)]
  )
  return(state_vectors)
}


#' find_adj_matrix
#'
#' Finds the adjacency matrix A for a given time window. For use within the
#' generate_adj_array function.
#'
#' @param state_vectors List of 2 state vectors x(t) and x(t+1) generated by
#'   generate_state_vectors
#' @param N Number of electrodes
#' @param t_window Length of one time window
#' @param nlambda Number of lambdas to be generated by glmnet cross-validation
#'
#' @return
#'
#' @examples
#' A is the 3D adjacency array in the following example, where k is the timewindow being calculated.
#' A[,,k] <- find_adj_matrix(state_vectors, N = 85, t_window = 250, nlambda = 16)
find_adj_matrix <- function(state_vectors, N, t_window, nlambda) {
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
  cv.ridge <- glmnet::cv.glmnet(H, b, alpha = 0, nfolds = 3, parallel = FALSE, nlambda = nlambda)
  lambdas <- rev(cv.ridge$lambda)
  
  test_lambda <- function(l, H, b) {
    ridge <- glmnet::glmnet(H, b, alpha = 0, lambda = l)
    N <- sqrt(dim(H)[2])
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
  stable_i <- FALSE
  
  while (!stable_i) {
    results <- test_lambda(lambdas[l], H = H, b = b)
    stable_i <- results$stable
    
    # if ((l+ncores-1) <= length(lambdas)) {
    #   results <- rave::lapply_async3(lambdas[l:(l+ncores-1)], test_lambda, H = H, b = b, .ncores = ncores)
    #   stable_i <- which(unname(unlist(lapply(results, function (x) x['stable']))))
    # } else {
    #   results <- rave::lapply_async3(lambdas[l:length(lambdas)], test_lambda, H = H, b = b, .ncores = ncores)
    #   stable_i <- which(unname(unlist(lapply(results, function (x) x['stable']))))
    #   break
    # }
    
    l <- l + 1
    
    if (l > length(lambdas)) {
      break
    }
  }
  
  if (!stable_i) {
    stop('No lambdas result in a stable adjacency matrix. Increase the number of lambdas, or (more likely) there is something wrong with your data.')
  }
  
  adj_matrix <- results$adj
  
  return(adj_matrix)
}

#' find_fragility 
#' 
#' Finds the fragility value for a single electrode (node) in a
#' single time window. For use within the generate_fragility_matrix function.
#'
#' @param node Integer specifying the node whose fragility is being calculated
#' @param A The adjacency matrix of the given time window
#' @param N Integer specifying number of electrodes
#'
#' @return
#'
#' @examples
#' 
find_fragility <- function(node, A_k, N, limit) {
  
  e_k <- vector(mode = 'numeric', length = N)
  e_k[node] <- 1
  
  argument <- t(e_k) %*% (solve(A_k - limit*diag(N))) # column perturbation
  # argument <- t(e_k) %*% t(solve(A_k - num*diag(N))) # row perturbation
  
  B <- rbind(Im(argument),Re(argument))
  
  perturb_mat <- (t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) %*% t(e_k) # column
  # perturb_mat <- e_k %*% t(t(B) %*% solve(B %*% t(B)) %*% c(0,-1)) # row
  
  norm(perturb_mat, type = '2')
}

#' draw_f_map_table
#'
#' Generates the parameters needed for the 3D brain viewer, fragility map,
#' most/least fragility table, and selected trials outputs.
#'
#' @param tnum An integer vector specifying which trials are selected. If
#'   multiple trials are selected, their fragility maps/values will be averaged.
#' @param adj_info The adj_info file generated by generate_adj_array.
#' @param f_path Pathname to subject's fragility files, with the trial number
#'   left empty. Will generally be something similar to
#'   'rave_data/data_dir/ProjectName/SubjectCode/rave/module_data/SubjectCode_f_info_trial_'
#' @param subject_code Character object specifying subject's code.
#' @param requested_electrodes Integer vector specifying which electrodes to
#'   display.
#' @param sort_fmap Specifies whether to sort map by electrode number or by
#'   fragility. Will be either "Electrode" or "Fragility".
#' @param check Check list generated by check_subject.
#' @param f_list_length Integer specifying how many values should be included on
#'   the most/least fragile list.
#'
#' @return
#' @export
#'
#' @examples
#' outputs <- draw_f_map_table(
#'     tnum = c(1,2,3), 
#'     adj_info = adj_info, 
#'     f_path = 'rave_data/data_dir/OnsetZone/PT01/rave/module_data/PT01_f_info_trial_', 
#'     subject_code = 'PT01', 
#'     requested_electrodes = c(1:24,26:36,42:43,46:54,56:70,72:95), 
#'     sort_fmap = 'Electrodes', 
#'     check = check, 
#'     f_list_length = 10
#' )
draw_f_map_table <- function(tnum, adj_info, f_path, subject_code, requested_electrodes, sort_fmap, check, f_list_length) {
  if (length(tnum) > 1) {
    f_plot <- list(
      norm = matrix(data = 0, nrow = dim(adj_info$A)[1], ncol = dim(adj_info$A)[3]),
      avg = vector(mode = 'numeric', length = dim(adj_info$A)[1]),
      trial = numeric()
    )
    for (i in seq_along(tnum)) {
      f_i <- readRDS(paste0(f_path,tnum[i]))
      f_plot[1:2] <- mapply(function(x,y) x + y, f_plot[1:2], f_i[2:3])
      f_plot$trial <- c(f_plot$trial,tnum[i])
    }
    f_plot[1:2] <- lapply(f_plot[1:2], function(x) x/length(tnum))
  } else {
    f <- readRDS(paste0(f_path,tnum))
    f_plot <- f[2:4]
  }
  
  f_plot$norm <- f_plot$norm[as.character(requested_electrodes),]
  f_plot$avg <- f_plot$avg[as.character(requested_electrodes)]
  brain_f <- data.frame("Subject"=subject_code,
                        "Electrode"=requested_electrodes,"Time"=0,
                        "Avg_Fragility"=f_plot$avg)
  
  elecsort <- sort(as.numeric(attr(f_plot$norm, "dimnames")[[1]]))
  fsort <- as.numeric(attr(sort(f_plot$avg), "names"))
  
  if (sort_fmap == 'Electrode') {
    elec_order <- elecsort
  } else if (sort_fmap == 'Fragility') {
    elec_order <- fsort
  }
  
  if (is.vector(f_plot$norm)){
    elec_order <- requested_electrodes
    x <- 1:length(f_plot$norm)
    m <- t(t(f_plot$norm))
  } else {
    f_plot$norm <- f_plot$norm[as.character(elec_order),]
    x <- 1:dim(f_plot$norm)[2]
    m <- t(f_plot$norm)
  }
  
  attr(m, 'xlab') = 'Time Window'
  attr(m, 'ylab') = 'Electrode'
  attr(m, 'zlab') = 'Fragility'
  
  if (check$elist) {
    y <- paste0(check$elec_list$Label[elec_order], '(', elec_order, ')')
    f_list <- paste0(check$elec_list$Label[fsort], '(', fsort, ')')
  } else {
    y <- elec_order
    f_list <- fsort
  }
  
  f_plot_params <- list(
    mat = m,
    x = x,
    y = y,
    zlim = c(0,1)
  )
  
  f_table_params <- data.frame(
    # Ranking = 1:f_list_length,
    Most.Fragile = rev(f_list)[1:f_list_length],
    Least.Fragile = f_list[1:f_list_length]
  )
  
  sel <- f_plot$trial
  
  output <- list(
    brain_f = brain_f,
    f_plot_params = f_plot_params,
    f_table_params = f_table_params,
    sel = sel
  )
}