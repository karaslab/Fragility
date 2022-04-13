# Put utility functions here, use `@export` to make function visible

#' @export
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