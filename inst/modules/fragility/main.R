# Main algorithm - rave_executes

require(Fragility)

# Initialize inputs
dev_Fragility(expose_functions = TRUE)

# mount_demo_subject()

rave::rave_prepare(
  subject = 'OnsetZone/PT01',
  electrodes =  c(1:24,26:36,42:43,46:54,56:70,72:95),
  epoch = 'PT01_sz',
  time_range = c(20,20),
  data_types = 'voltage',
  reference = 'car'
)

init_module('fragility', debug = TRUE)

# rave::rave_prepare(
#   subject = 'OnsetZone/KAA',
#   electrodes =  c(1:43,45,46,61:72,89:116,145:158,196:209),
#   epoch = 'KAA_sz',
#   time_range = c(20,20),
#   data_types = 'voltage',
#   reference = 'car'
# )

# >>>>>>>>>>>> Start ------------- [DO NOT EDIT THIS LINE] ---------------------
######' @auto=TRUE

# only use electrodes and trials requested
requested_electrodes = dipsaus::parse_svec(text_electrode)

trial = module_tools$get_meta('trials')

# trial for calculating adj matrix
tnum_adj <- trial$Trial[trial$Condition %in% adj_conditions]

# trial(s) for display on fragility map
tnum <- trial$Trial[trial$Condition %in% requested_conditions]

# where subject's pt_info, adj_info, and f_info RDS files are saved
subject_dir <- module_tools$get_subject_dirs()$module_data_dir
if (!file.exists(subject_dir)){
  dir.create(subject_dir)
}
subject_code <- subject$subject_code

srate <- module_tools$get_sample_rate(original = TRUE)

# check if subject has pt_info, adj_info, and f_info files
check <- check_subject(subject_code,subject_dir,trial$Trial)

if (check$pt) {
  pt_info_all <- readRDS(paste0(subject_dir,'/',subject_code,'_pt_info'))
  pt_info <- list(
    v = pt_info_all$v[tnum,,as.character(requested_electrodes)],
    trial = tnum
  )
} else {
  print("pt has not previously been loaded. click load patient button")
}

if (check$adj[tnum_adj]) {
  adj_info <- readRDS(paste0(subject_dir,'/',subject_code,'_adj_info_trial_',tnum_adj))
} else {
  print('adj matrix for this trial doesnt exist yet. click generate adjacency matrix button')
}

S <- dim(pt_info_all$v)[2] # S is total number of timepoints

if(S %% requested_tstep != 0) {
  # truncate S to greatest number evenly divisible by timestep
  S <- trunc(S/requested_tstep) * requested_tstep
}
J <- S/requested_tstep - (requested_twindow/requested_tstep) + 1

t_estimate <- 1.5*J
t_estimate_hrs_min <- paste0(t_estimate%/%60, ' hrs, ', round(t_estimate%%60, digits = 1), ' min')
message_board <- paste0('Estimated time: ', t_estimate_hrs_min)

if (all(check$f[tnum])) {
  if (length(tnum) > 1) {
    f_norm_plot <- matrix(data = 0, nrow = dim(adj_info$A)[1], ncol = dim(adj_info$A)[3])
    f_avg_plot <- vector(mode = 'numeric', length = dim(adj_info$A)[1])
    for (i in tnum) {
      f_info <- readRDS(paste0(subject_dir,'/',subject_code,'_f_info_trial_',tnum[i]))
      f_norm_plot <- f_info$norm + f_norm_plot
      f_avg_plot <- f_info$avg + f_avg_plot
    }
    f_norm_plot <- f_norm_plot/length(tnum)
    f_avg_plot <- f_avg_plot/length(tnum)
  } else {
    f_info <- readRDS(paste0(subject_dir,'/',subject_code,'_f_info_trial_',tnum))
    f_norm_plot <- f_info$norm
    f_avg_plot <- f_info$avg
  }
  f_norm_plot <- f_norm_plot[as.character(requested_electrodes),]
  f_avg_plot <- f_avg_plot[as.character(requested_electrodes)]
  
  if (sort_fmap == 'Electrode') {
    elecsort <- sort(as.numeric(attr(f_norm_plot, "dimnames")[[1]]))
    f_norm_plot <- f_norm_plot[as.character(elecsort),]
  } else if (sort_fmap == 'Fragility') {
    elecsort <- as.numeric(attr(sort(f_avg_plot), "names"))
    f_norm_plot <- f_norm_plot[as.character(elecsort),]
  }
  
  f_plot_params <- list(
    mat = f_norm_plot,
    x = 1:dim(f_norm_plot)[2],
    y = 1:dim(f_norm_plot)[1],
    zlim = c(0,1)
  )
} else {
  print('fragility map for some trials has not been generated yet. click generate fragility matrix button')
}

# <<<<<<<<<<<< End ----------------- [DO NOT EDIT THIS LINE] -------------------

# Debug

Fragility::dev_Fragility(expose_functions = TRUE)

# Debug - offline:
main = Fragility:::debug_module('fragility')
ret = main()
result = ret$results

result$get_value('preload_info')


# Debug - online:
Fragility::dev_Fragility(expose_functions = TRUE)
# mount_demo_subject()
rave::rave_prepare(
  subject = 'OnsetZone/PT01',
  electrodes =  c(1:24,26:36,42:43,46:54,56:70,72:95),
  epoch = 'PT01_sz',
  time_range = c(20,20),
  data_types = 'voltage',
  reference = 'car'
)
view_layout('fragility')

# Production - Deploy as RAVE module
# Always Ctrl/cmd+shift+B first if you want to deploy it online
rm(list = ls(all.names = TRUE)); rstudioapi::restartSession()
module = rave::get_module(package = 'Fragility', module_id = 'fragility')
rave::init_app(module)
