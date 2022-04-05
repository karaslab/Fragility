# Put utility functions here, use `@export` to make function visible

#' @export
check_subject <- function(subject_code, subject_dir, trials) {
  
  pt <- file.exists(paste0(subject_dir,'/',subject_code,'_pt_info'))
  
  adj <- vector(length = length(trials))
  f <- vector(length = length(trials))
  
  for (i in trials) {
    adj[i] <- file.exists(paste0(subject_dir,'/',subject_code,'_adj_info_trial_',i))
    f[i] <- file.exists(paste0(subject_dir,'/',subject_code,'_f_info_trial_',i))
  }
  
  check <- list(
    pt = pt,
    adj = adj,
    f = f
  )
}