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

# initalize pt_info if available
if (local_data$check$pt && is.null(isolate(local_data$pt_info))) {
  print('Loading pt_info first time')
  showNotification('Loading existing pre-processed patient...', id = 'pt_loading')
  # load it in if the file exists but hasn't been loaded in yet
  local_data$pt_info <- readRDS(paste0(module_data,'/',subject_code,'_pt_info'))
  removeNotification('pt_loading')
  if (!all(as.character(preload_info$electrodes) %in% attr(local_data$pt_info$v, "dimnames")$Electrode)){
    stop('Not all loaded electrodes are present in saved pt_info! Please reload 
         the patient selecting only electrodes that were previously saved. 
         Alternatively, re-process the patient to include all loaded electrodes.')
  }
}

# initialize adj_info if available
if (any(local_data$check$adj) & (is.null(isolate(local_data$adj_info)) | new_subject)) {
  print('Loading adj_info first time')
  showNotification('Loading existing adjacency array...', id = 'adj_loading')
  tnum_adj <- which(local_data$check$adj)[1]
  local_data$adj_info <- readRDS(paste0(module_data,'/',subject_code,'_adj_info_trial_',tnum_adj))
  local_data$selected$adj <- local_data$adj_info$trial
  removeNotification('adj_loading')
  new_subject <- FALSE
} else {
  # if already initialized, update tnum_adj to reflect user selection
  tnum_adj <- trial$Trial[trial$Condition %in% adj_conditions]
}

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
  electrodes =  c(1:24,26:36,42:43,46:54,56:70,72:95),
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

view_layout('fragility')

# Production - Deploy as RAVE module
# Always Ctrl/cmd+shift+B first if you want to deploy it online
rm(list = ls(all.names = TRUE)); rstudioapi::restartSession()
module = rave::get_module(package = 'Fragility', module_id = 'fragility')
rave::init_app(module)

result = module()
