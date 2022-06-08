# Main algorithm - rave_executes

require(Fragility)

# Initialize inputs
dev_Fragility(expose_functions = TRUE)

# mount_demo_subject()

# mount PT01 as demo subject
rave::rave_prepare(
  subject = 'OnsetZone/PT01',
  electrodes =  c(1:24,26:36,42:43,46:54,56:70,72:95),
  epoch = 'PT01_sz',
  time_range = c(20,20),
  data_types = 'voltage',
  reference = 'car'
)

# initiates module in R environment
init_module('fragility', debug = TRUE)

# >>>>>>>>>>>> Start ------------- [DO NOT EDIT THIS LINE] ---------------------
######' @auto=TRUE

print('executing main')

# initialize and load pt_info if available
local_data$elec_present <- TRUE
if (local_data$check$pt & (is.null(isolate(local_data$pt_info)) | new_subject)) {
  print('Loading pt_info first time')
  showNotification('Loading existing pre-processed patient...', id = 'pt_loading')
  
  # load it in if the file exists but hasn't been loaded in yet
  local_data$pt_info <- readRDS(paste0(module_data,'/',subject_code,'_pt_info'))
  removeNotification('pt_loading')
  
  # check if saved pt_info electrodes match the loaded electrodes
  elec_check <- as.character(preload_info$electrodes) %in% attr(local_data$pt_info$v, "dimnames")$Electrode
  if (!all(elec_check)){
    missing_i <- which(!elec_check)
    
    local_data$elec_present <- FALSE
    
    shiny::showNotification(paste0('Not all loaded electrodes are present in saved pt_info! 
                Electrodes missing from pt_info: ', 
                dipsaus::deparse_svec(preload_info$electrodes[missing_i])), duration = NULL, type = 'error')
    shiny::showNotification('Please reload the patient selecting only electrodes that were previously saved. 
                            Alternatively, re-process the patient to include all loaded electrodes.', duration = NULL, type = 'error')
  }
  
  if (!any(local_data$check$adj)) {
    new_subject <- FALSE
  }
}

# initialize adj_info if available
if (any(local_data$check$adj) & (is.null(isolate(local_data$adj_info)) | new_subject)) {
  print('Loading adj_info first time')
  showNotification('Loading existing adjacency array...', id = 'adj_loading')
  
  # automatically load first trial as default
  tnum_adj <- which(local_data$check$adj)[1]
  local_data$adj_info <- readRDS(paste0(module_data,'/',subject_code,'_adj_info_trial_',tnum_adj))
  local_data$selected$adj <- local_data$adj_info$trial
  removeNotification('adj_loading')
  new_subject <- FALSE
} else {
  # if already initialized, update tnum_adj to reflect user selection
  tnum_adj <- trial$Trial[trial$Condition %in% adj_conditions]
}

# check which files are available every recalculate
local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)

# update available trials every recalculate
possible <- list(
  pt = local_data$check$pt,
  adj = trial$Trial[local_data$check$adj],
  f = trial$Trial[local_data$check$f]
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
  electrodes =  c(1:4,7:24,26:36,42:43,46:54,56:69,72:95),
  epoch = 'PT01_sz',
  time_range = c(20,20),
  data_types = 'voltage',
  reference = 'car'
)

rave::rave_prepare(
  subject = 'OnsetZone/KAA',
  electrodes =  c(1:43,45:80,82:116,129:164,166:184,186:223,226:235,237:244),
  epoch = 'KAA_sz',
  time_range = c(20,20),
  data_types = 'voltage',
  reference = 'car'
)

rave::rave_prepare(
  subject = 'OnsetZone/YDS_seizure',
  electrodes =  c(1:19,21:35,42:46,48:219),
  epoch = 'YDS_seizure',
  time_range = c(20,20),
  data_types = 'voltage',
  reference = 'car'
)

rave::rave_prepare(
  subject = 'OnsetZone/YDO_seizure',
  electrodes =  c(1:19,21:37,44:241),
  epoch = 'YDO_seizure',
  time_range = c(20,20),
  data_types = 'voltage',
  reference = 'car'
)

rave::rave_prepare(
  subject = 'OnsetZone/PT026',
  electrodes =  c(1:4,7:12,15:23,25:33,47:63,65:66,69:71,73:110),
  epoch = 'PT_026_sz',
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
