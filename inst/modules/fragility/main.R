# Main algorithm - rave_executes

require(Fragility)

# Initialize inputs
dev_Fragility(expose_functions = TRUE)

mount_demo_subject()

init_module('fragility', debug = TRUE)

rave::rave_prepare(
  subject = 'KAA',
  electrodes =  c(1:43,45:46,61:72,89:116,145:158,196:209),
  epoch = 'KAA_sz',
  time_range = c(20,20),
  data_types = 'voltage',
  reference = 'car'
)

# >>>>>>>>>>>> Start ------------- [DO NOT EDIT THIS LINE] ---------------------
######' @auto=TRUE

requested_electrodes = dipsaus::parse_svec(text_electrode)

time_points = preload_info$time_points
frequencies = preload_info$frequencies

trial = module_tools$get_meta('trials')

# use only the electrode(s) and trial type(s) requested

tnum <- trial$Trial[trial$Condition %in% requested_conditions]

# load the power
p <- module_tools$get_power()

# subset based on UI choices
p.sub <- p$subset(
  Trial ~ trial$Condition %in% requested_conditions,
  Electrode ~ Electrode %in% requested_electrodes,
  Frequency ~ Frequency %within% requested_frequencies,
  # just return an array, rather than a Tensor object
  data_only = TRUE,
  #don't drop length=1 margins
  drop = FALSE
)

# calculate baseline-corrected data across time(along_dim = 3), per frequency, per trial, per electrode; (unit_dims = 1,2,4)
bsl <- dipsaus::baseline_array(
  x=p.sub, along_dim = 3, unit_dims = c(1,2,4), method = 'percentage',
  baseline_indexpoints = which(
    time_points >= min(requested_baseline) & 
      time_points  <= max(requested_baseline))
)

pt_info <- load_fragility_patient(
  subject_code = 'KAA',
  block = '4',
  elec = requested_electrodes
)

adj_info <- generate_adj_array(
  t_window = 250, 
  t_step = 125, 
  v = pt_info$v, 
  ncores = 8
)

f_info <- generate_fragility_matrix(
  N = pt_info$N,
  J = adj_info$J,
  A = adj_info$A,
  elec = requested_electrodes
)

f_plot_params <- list(
  mat = f_info$norm,
  x = 1:J,
  y = 1:N,
  zlim = c(0,1)
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
mount_demo_subject()
view_layout('fragility')

# Production - Deploy as RAVE module
# Always Ctrl/cmd+shift+B first if you want to deploy it online
rm(list = ls(all.names = TRUE)); rstudioapi::restartSession()
module = rave::get_module(package = 'Fragility', module_id = 'fragility')
rave::init_app(module)
