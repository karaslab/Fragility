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

print('executing main')
# print(text_electrode)
# print(recording_unit)
# print(adj_conditions)
# print(requested_conditions)
# print(requested_twindow)
# print(requested_tstep)
# print(requested_nlambda)
# print(requested_ncores)
# print(future_maxsize)
# print(sort_fmap)
# only use electrodes and trials requested
requested_electrodes = dipsaus::parse_svec(text_electrode)

# trial for calculating adj and fragility matrix
tnum_adj <- trial$Trial[trial$Condition %in% adj_conditions]

# trial(s) for display on fragility map
tnum <- trial$Trial[trial$Condition %in% requested_conditions]

# where subject's pt_info, adj_info, and f_info RDS files are saved
# subject_dir <- module_tools$get_subject_dirs()$module_data_dir
# if (!file.exists(subject_dir)){
#   dir.create(subject_dir)
# }
# subject_code <- subject$subject_code
# 
# srate <- module_tools$get_sample_rate(original = TRUE)

# check if subject has pt_info, adj_info, and f_info files
# check <- check_subject(subject_code,subject_dir,trial$Trial)
# print(input$requested_conditions %in% module_tools$get_meta('trials')$Condition[check$f])
# updateSelectInput(session = session, inputId = 'requested_conditions',
#                   choices = module_tools$get_meta('trials')$Condition[check$f],
#                   )
# print(str(input$requested_conditions))
# print(str(module_tools$get_meta('trials')$Condition[check$f]))
#                   # selected = input$requested_conditions)


if (check$pt) {
  pt_info_all <- readRDS(paste0(subject_dir,'/',subject_code,'_pt_info'))
  # pt_info <- list(
  #   v = pt_info_all$v[tnum,,as.character(requested_electrodes)],
  #   trial = tnum
  # )
  
  # Show estimated adj array calculation time
  if (requested_tstep != "") {
    S <- dim(pt_info_all$v)[2] # S is total number of timepoints
    
    tstep <- as.numeric(requested_tstep)
    
    if(S %% tstep != 0) {
      # truncate S to greatest number evenly divisible by timestep
      S <- trunc(S/tstep) * tstep
    }
    J <- S/tstep - (requested_twindow/tstep) + 1
    
    local_data$J <- J
    
    t_estimate <- 1.5*J
    t_estimate_hrs_min <- paste0(t_estimate%/%60, ' hours, ', round(t_estimate%%60, digits = 1), ' minutes')
    local_data$estimate <- paste0('Estimated time: ', t_estimate_hrs_min)
  }
  
} else {
  print("pt has not previously been loaded. click load patient button")
}

if (check$adj[tnum_adj]) {
  adj_info <- readRDS(paste0(subject_dir,'/',subject_code,'_adj_info_trial_',tnum_adj))
} else {
  print('adj matrix for this trial doesnt exist yet. click generate adjacency matrix button')
}

# fragility map stuff
# do.call(req, split(check$f, seq_along(check$f)))
if (length(tnum) > 1) {
  f_plot <- list(
    norm = matrix(data = 0, nrow = dim(adj_info$A)[1], ncol = dim(adj_info$A)[3]),
    avg = vector(mode = 'numeric', length = dim(adj_info$A)[1]),
    trial = numeric()
  )
  for (i in seq_along(tnum)) {
    f_info <- readRDS(paste0(subject_dir,'/',subject_code,'_f_info_trial_',tnum[i]))
    f_plot[1:2] <- mapply(function(x,y) x + y, f_plot[1:2], f_info[2:3])
    # f_plot$norm <- f_plot$norm + f_info$norm
    # f_plot$avg <- f_plot$avg + f_info$avg
    f_plot$trial <- c(f_plot$trial,tnum[i])
  }
  f_plot[1:2] <- lapply(f_plot[1:2], function(x) x/length(tnum))
} else {
  f_info <- readRDS(paste0(subject_dir,'/',subject_code,'_f_info_trial_',tnum))
  f_plot <- f_info[2:4]
}

f_plot$norm <- f_plot$norm[as.character(requested_electrodes),]
f_plot$avg <- f_plot$avg[as.character(requested_electrodes)]
local_data$brain_f <- data.frame("Subject"=subject_code,
                                 "Electrode"=requested_electrodes,"Time"=0,
                                 "Avg_Fragility"=f_plot$avg)

if (sort_fmap == 'Electrode') {
  elecsort <- sort(as.numeric(attr(f_plot$norm, "dimnames")[[1]]))
  f_plot$norm <- f_plot$norm[as.character(elecsort),]
} else if (sort_fmap == 'Fragility') {
  elecsort <- as.numeric(attr(sort(f_plot$avg), "names"))
  f_plot$norm <- f_plot$norm[as.character(elecsort),]
}

f_plot_params <- list(
  mat = f_plot$norm,
  x = 1:dim(f_plot$norm)[2],
  y = 1:dim(f_plot$norm)[1],
  zlim = c(0,1)
)

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
view_layout('fragility')

# Production - Deploy as RAVE module
# Always Ctrl/cmd+shift+B first if you want to deploy it online
rm(list = ls(all.names = TRUE)); rstudioapi::restartSession()
module = rave::get_module(package = 'Fragility', module_id = 'fragility')
rave::init_app(module)
