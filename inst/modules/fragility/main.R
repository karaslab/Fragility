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

rave::rave_prepare(
  subject = 'OnsetZone/KAA',
  electrodes =  c(1:116,129:244),
  epoch = 'KAA_sz',
  time_range = c(20,20),
  data_types = 'voltage',
  reference = 'car'
)

# >>>>>>>>>>>> Start ------------- [DO NOT EDIT THIS LINE] ---------------------
######' @auto=TRUE

print('executing main')
req(adj_conditions, requested_conditions)

# only use electrodes and trials requested
requested_electrodes = dipsaus::parse_svec(text_electrode)

# trial for calculating adj and fragility matrix
tnum_adj <- trial$Trial[trial$Condition %in% adj_conditions]

# temporary fixes for RAVE initialization errors (inputs being NULL)

# if (is.null(adj_conditions)) {
#   tnum_adj <- 1
# }
# 
# if(is.null(requested_conditions)) {
#   tnum <- 1
# }
# 
# if(is.null(requested_tstep)) {
#   requested_tstep <- ""
# }

# trial(s) for display on fragility map
tnum <- trial$Trial[trial$Condition %in% requested_conditions]

# # where subject's pt_info, adj_info, and f_info RDS files are saved
# subject_dir <- module_tools$get_subject_dirs()
# module_data <- subject_dir$module_data_dir
# if (!file.exists(subject_dir$module_data_dir)){
#   dir.create(subject_dir$module_data_dir)
# }
# subject_code <- subject$subject_code
# 
# srate <- module_tools$get_sample_rate(original = TRUE)

# local_data$check if subject has pt_info, adj_info, and f_info files
# local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
local_data$check <- check
# print(str(local_data$check))

if (local_data$check$pt) {
  pt_info <- readRDS(paste0(module_data,'/',subject_code,'_pt_info'))
  
  if (local_data$check$adj[tnum_adj]) {
    adj_info <- readRDS(paste0(module_data,'/',subject_code,'_adj_info_trial_',tnum_adj))
  } else {
    adj_info <- list(trial = "None")
    print('adj matrix for this trial doesnt exist yet. click generate adjacency matrix button')
    showNotification('No valid adjacency arrays detected. Please choose a trial and click the "Generate Adjacency Matrix" button.', duration = 10)
  }
} else {
  adj_info <- list(trial = "None")
  print("pt has not previously been loaded. click load patient button")
  showNotification('Patient has not been processed yet. Please click the "Pre-process Patient" button under "Load Patient".', duration = 10)
}

# fragility map stuff
# do.call(req, split(local_data$check$f, seq_along(local_data$check$f)))
if (any(local_data$check$f)){
  if (length(tnum) > 1) {
    f_plot <- list(
      norm = matrix(data = 0, nrow = dim(adj_info$A)[1], ncol = dim(adj_info$A)[3]),
      avg = vector(mode = 'numeric', length = dim(adj_info$A)[1]),
      trial = numeric()
    )
    for (i in seq_along(tnum)) {
      f_info <- readRDS(paste0(module_data,'/',subject_code,'_f_info_trial_',tnum[i]))
      f_plot[1:2] <- mapply(function(x,y) x + y, f_plot[1:2], f_info[2:3])
      # f_plot$norm <- f_plot$norm + f_info$norm
      # f_plot$avg <- f_plot$avg + f_info$avg
      f_plot$trial <- c(f_plot$trial,tnum[i])
    }
    f_plot[1:2] <- lapply(f_plot[1:2], function(x) x/length(tnum))
  } else {
    f_info <- readRDS(paste0(module_data,'/',subject_code,'_f_info_trial_',tnum))
    f_plot <- f_info[2:4]
  }
  
  f_plot$norm <- f_plot$norm[as.character(requested_electrodes),]
  f_plot$avg <- f_plot$avg[as.character(requested_electrodes)]
  local_data$brain_f <- data.frame("Subject"=subject_code,
                                   "Electrode"=requested_electrodes,"Time"=0,
                                   "Avg_Fragility"=f_plot$avg)
  
  elecsort <- sort(as.numeric(attr(f_plot$norm, "dimnames")[[1]]))
  fsort <- as.numeric(attr(sort(f_plot$avg), "names"))
  
  if (sort_fmap == 'Electrode') {
    elec_order <- elecsort
  } else if (sort_fmap == 'Fragility') {
    elec_order <- fsort
  }
  
  f_plot$norm <- f_plot$norm[as.character(elec_order),]
  
  m <-  t(f_plot$norm)
  attr(m, 'xlab') = 'Time'
  attr(m, 'ylab') = 'Electrode'
  attr(m, 'zlab') = 'Fragility'

  if (local_data$check$elist) {
    y <- paste0(local_data$check$elec_list$Label[elec_order], '(', elec_order, ')')
    f_list <- paste0(local_data$check$elec_list$Label[fsort], '(', fsort, ')')
  } else {
    y <- elec_order
    f_list <- fsort
  }
  
  f_plot_params <- list(
    mat = m,
    x = 1:dim(f_plot$norm)[2],
    y = y,
    zlim = c(0,1)
  )
  
  f_table_params <- data.frame(
    # Ranking = 1:f_list_length,
    Most.Fragile = rev(f_list)[1:f_list_length],
    Least.Fragile = f_list[1:f_list_length]
  )
  
} else {
  f_plot <- list(trial = 'None')
}

selected <- list(
  adj = adj_info$trial,
  f = f_plot$trial
)


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

rave::rave_prepare(
  subject = 'OnsetZone/KAA',
  electrodes =  c(1:116,129:244),
  epoch = 'KAA_sz',
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

result = module()
