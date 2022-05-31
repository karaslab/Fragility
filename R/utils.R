# Put utility functions here, use `@export` to make function visible

#' check_subject
#'
#' @param subject_code Character specifying subject's code
#' @param subject_dir List containing sub-directories of subject's folder in
#'   rave_data, from module_tools$get_subject_dirs()
#' @param trials Integer vector with a number for each trial
#'
#' @export
#' @examples
#' check <- check_subject('KAA', module_tools$get_subject_dirs(), 1:8)
check_subject <- function(subject_code, subject_dir, trials) {
  
  module_data <- subject_dir$module_data_dir
  
  pt <- file.exists(paste0(module_data,'/',subject_code,'_pt_info'))
  
  adj <- vector(length = length(trials))
  f <- vector(length = length(trials))
  
  for (i in trials) {
    adj[i] <- file.exists(paste0(module_data,'/',subject_code,'_adj_info_trial_',i))
    f[i] <- file.exists(paste0(module_data,'/',subject_code,'_f_info_trial_',i))
  }
  
  if (file.exists(paste0(subject_dir$meta_dir,'/electrodes.csv'))){
    elec_list <- read.csv(paste0(subject_dir$meta_dir,'/electrodes.csv'))
    elist <- !all(elec_list$Label == 'NoLabel')
  } else {
    elec_list <- NULL
    elist <- FALSE
  }
  
  check <- list(
    pt = pt,
    adj = adj,
    f = f,
    elec_list = elec_list,
    elist = elist
  )
}

# estimate_time <- function(pt_info, tstep, twindow) {
#   S <- dim(pt_info$v)[2] # S is total number of timepoints
#   N <- dim(pt_info$v)[3]
#   
#   if(S %% tstep != 0) {
#     # truncate S to greatest number evenly divisible by timestep
#     S <- trunc(S/tstep) * tstep
#   }
#   J <- S/tstep - (twindow/tstep) + 1
#   
#   t_estimate <- 1.5*J
#   
#   
#   # est_time <- paste0('Estimated time: ', t_estimate%/%60, ' hours, ', round(t_estimate%%60, digits = 1), ' minutes')
#   # warning this H can be big
#   # Hsize <- object.size(matrix(0, nrow = N*(twindow-1), ncol = N^2)) + (200 * 1024^2)
#   
#   time_estimate <- list(
#     time = est_time,
#     # Hsize = Hsize,
#     J = J
#   )
# }

