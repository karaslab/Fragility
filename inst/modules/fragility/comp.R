# File defining module inputs, outputs

# ----------------------------------- Debug ------------------------------------
require(shiny)
require(rave)
require(Fragility)

env = dev_Fragility(TRUE)

#' Load subject for debugging
#' Make sure this is executed before developing the module to make sure
#' at least one subject is loaded
mount_demo_subject()


# >>>>>>>>>>>> Start ------------- [DO NOT EDIT THIS LINE] ---------------------


#  ----------------------  Initializing Global variables -----------------------

#' define_initialization is executed every time when:
#'   *** a subject is loaded
#'
#' You can use this function to do:
#'   1. Check if data is complete
#'   2. Define some "global" variables
#'
#' Check package vignette to see the explanation

define_initialization({
  # Enter code to handle data when a new subject is loaded
  print('initializing')
  
  # load voltage data
  rave_checks('voltage referenced')
  volt <- module_tools$get_voltage()
  local_data$v <- volt$get_data()
  
  # get necessary directories and subject data for file storage
  subject_dir <- module_tools$get_subject_dirs()
  module_data <- subject_dir$module_data_dir
  
  # if module_data directory doesn't exist yet, create it
  if (!file.exists(module_data)){
    dir.create(module_data)
  }

  # get various subject metadata
  subject_code <- subject$subject_code
  trial <- module_tools$get_meta('trials')
  srate <- module_tools$get_sample_rate(original = TRUE)
  
  # check if subject has pt_info, adj_info, and f_info files
  local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
  
  new_subject <- TRUE
})

load_scripts(
  'inst/modules/fragility/event_handlers.R', 
  asis = TRUE
)


#  ---------------------------------  Inputs -----------------------------------
#' Define inputs
#'
#' Use function `define_input` to define inputs.
#' Make sure rave toolbox (dev_[your_package_name]) is loaded
#'
#' @Usage: define_input(definition, init_args, init_expr)
#'
#' @param definition defines input types, for example:
#'     textInput(inputId = 'text_electrode', label = 'Electrode')
#'   defines a character variable `text_electrode` as one of the inputs
#'
#' Here are some possible types
#'   1. textInput is an input for characters
#'   2. numericInput is for numbers
#'   3. sliderInput is for number or number ranges
#'   4. selectInput is to provide choices
#'   5. checkboxInput is for Yes/No options
#'   6. compoundInput is for multiple group inputs
#'   7. customizedUI is for advanced UI controls
#'
#'   Most of basic UI widgets can be found at:
#'       https://shiny.rstudio.com/gallery/widget-gallery.html
#'
#'   The easiest way to look for usage is using `help` function:
#'   help('textInput'), or ??textInput
#'
#'
#' @param init_args,init_expr define the action when a subject is loaded
#'
#' Use the definition
#'
#'     textInput(inputId = 'text_electrode', label = 'Electrode')
#'
#' as an example. You might want label to update according to electrodes loaded.
#' In this case, you can assign
#'     init_args = 'label'
#' indicating the label of this input will update according to actual data.
#'
#' You may also change multiple arguments. The following code is an extended example

# Step 1: Process Patient----

define_input(
  definition = selectInput(inputId = 'recording_unit', label='Units of Recording', choices = character(0)),
  
  init_args = c('choices','selected'),
  
  init_expr = {
    choices = c('mV','uV','nV')
    selected = 'uV'
  }
)

define_input(
  definition = checkboxInput(inputId = 'half_hz', label='Halve Sampling Rate?')
)

define_input(
  actionButton(inputId = 'process_pt', label='Pre-process Patient')
)

# Step 2: Generate Model----

define_input(
  sliderInput(inputId = 'requested_twindow', label='Time Window Size (ms)', min=100, max=1000, value=100, step=10, ticks = FALSE)
)

define_input(
  selectInput(inputId = 'requested_tstep', label='Time Step (ms)', choices = character(0)),
  
  init_args = c('choices', 'selected'),
  
  init_expr = {
    # set timestep choices to either size of timewindow or half of timewindow
    choices = c(input$requested_twindow, input$requested_twindow/2)
    
    # default timestep is equal to timewindow
    selected = input$requested_twindow
  }
)

define_input(
  numericInput(inputId = 'requested_nlambda', label='Number of Lambdas', min=1, max=100, value=16, step=1)
)

define_input(
  numericInput(inputId = 'requested_ncores', label='Number of cores', min=1, max=4, value=4, step=1),
  
  init_args = c('max', 'value'),
  
  init_expr = {
    # default number of cores to half of maximum
    max = as.numeric(rave::rave_options('max_worker'))
    value = (max + 1)/2
  }
)

define_input(
  selectInput(inputId = 'adj_conditions', choices = character(0), label='Seizure Trial for Matrix Generation'),
  
  init_args = c('choices', 'selected'),
  
  init_expr = {
    choices = unique(module_tools$get_meta('trials')$Condition)
    selected = module_tools$get_meta('trials')$Condition[1]
  }
)

define_input(
  actionButton(inputId = 'gen_adj', label='Generate Adjacency Matrix for this Trial')
)

define_input(
  actionButton(inputId = 'gen_f', label='Generate Fragility Matrix for this Trial')
)

# Step 3: View Fragility----

define_input(
  selectInput(inputId = 'requested_conditions', choices = character(0), multiple = TRUE, label='Seizure Trial(s) for Fragility Display'),
  
  init_args = c('label','choices', 'selected'),
  
  init_expr = {
    if (any(local_data$check$f)){
      # allow user to choose between trials that have available f_info files
      label = 'Seizure Trial(s) for Fragility Map Display'
      choices = module_tools$get_meta('trials')$Condition[local_data$check$f]
      selected = module_tools$get_meta('trials')$Condition[local_data$check$f][1]
    } else {
      # if there are no available f_info files
      label = 'No valid fragility matrices detected! Please generate one above.'
      choices = character(0)
      selected = NULL
    }
  }
)

define_input(
  definition = textInput(inputId = 'text_electrode', label = 'Electrode'),

  # label will tell users which electrodes are loaded
  # value will be the first electrode

  init_args = c('label', 'value'),
  init_expr = {

    # check ?rave_prepare for what kind of data are loaded
    loaded_electrodes = preload_info$electrodes

    # Generate text for loaded electrodes
    text = dipsaus::deparse_svec(loaded_electrodes)

    # Update
    label = paste0('Select Fragility Electrode (Loaded: ', text, ')')
    value = text
  }
)

define_input(
  selectInput(inputId = 'sort_fmap', choices = c('Electrode (descending)', 'Electrode (ascending)','Fragility (descending)', 'Fragility (ascending)'), selected = 'Electrode (descending)', label='Sort Fragility Map By:')
)

define_input(
  numericInput(inputId = 'f_list_length', label='List how many for most/least fragile?', min=1, max=1, value=1),
  init_args = c('max','value'),
  
  init_expr = {
    # default length of list is half of number of loaded electrodes
    max = floor(length(preload_info$electrodes)/2)
    value = min(c(floor(length(preload_info$electrodes)/8), 10))
  }
)

define_input(
  sliderInput(inputId = 'sz_onset', label='Seizure Onset Marker', min=1, max=1, value=1, step=1),
  
  init_args = c('min','max','value'),
  
  init_expr = {
    min = preload_info$time_points[1]
    max = preload_info$time_points[length(preload_info$time_points)]
    value = 0
  }
)

define_input(
  sliderInput(inputId = 'exponentiate', label='Exponentiate 3D Viewer Fragility (for higher contrast)', min=1, max=5, value=1, step=2, ticks = FALSE)
)

define_input(
  checkboxInput(inputId = 'exp_fmap', label = 'Exponentiate Fragility Map too?')
)

define_input(
  checkboxInput(inputId = 'auto_calc', label='Auto-recalculate?', value = TRUE)
)

define_input(
  actionButton(inputId = 'draw_f_map', label='Calculate Fragility!'),
  
  init_args = c('label'),
  
  init_expr = {
    
    # enable calculate fragility button only if fragility matrix files exist
    if (any(local_data$check$f)){
      label = 'Calculate Fragility!'
      shinyjs::enable('draw_f_map')
    } else {
      label = 'No valid fragility matrices detected!'
      shinyjs::disable('draw_f_map')
    }
  }
)

# View Original Raw EEG Data----

define_input(
  selectInput(inputId = 'v_conditions', choices = character(0), label='Seizure Trial for EEG Trace Viewer'),
  
  init_args = c('choices', 'selected'),
  
  init_expr = {
    choices = unique(module_tools$get_meta('trials')$Condition)
    selected = module_tools$get_meta('trials')$Condition[1]
  }
)

define_input(
  definition = textInput(inputId = 'v_electrode', label = 'Electrode'),
  
  init_args = c('label', 'value'),
  init_expr = {
    
    loaded_electrodes = preload_info$electrodes
    
    text = dipsaus::deparse_svec(loaded_electrodes)
    
    # Update
    label = paste0('Select EEG Trace Electrodes (Loaded: ', text, ')')
    value = text
  }
)

define_input(
  checkboxInput(inputId = 'v_sync', label = 'Sync electrodes with Fragility Map?', value = TRUE)
)

define_input(
  actionButton(inputId = 'load_v_traces', label='Load EEG Voltage Traces')
)

# Re-check Files----

define_input(
  actionButton(inputId = 'refresh_btn', label='Refresh')
)

# Experimental----
define_input(
  checkboxInput(inputId = 'experimental', label = 'Enable experimental features?')
)

define_input(
  numericInput(inputId = 'limreal', label = 'Real part of limit for fragility eigenvalue', value = 0)
)

define_input(
  numericInput(inputId = 'limimag', label = 'Imaginary part of limit for fragility eigenvalue', value = 1)
)

# the input_layout list is used by rave to determine order and grouping of layouts
input_layout <- list(
  '[-]Step 1: Process Patient' = list(
    'recording_unit',
    'half_hz',
    'process_pt'
  ),
  '[-]Step 2: Generate Model' = list(
    'requested_twindow',
    'requested_tstep',
    'requested_nlambda',
    'requested_ncores',
    'adj_conditions',
    'gen_adj',
    'gen_f'
  ),
  '[-]Step 3: View Fragility' = list(
    'requested_conditions',
    'text_electrode',
    'sort_fmap',
    'f_list_length',
    'sz_onset',
    'exponentiate',
    'exp_fmap',
    'auto_calc',
    'draw_f_map'
  ),
  '[-]View Original Raw EEG Data' = list(
    'v_conditions',
    'v_electrode',
    'v_sync',
    'load_v_traces'
  ),
  '[-]Re-check Files' = list(
    'refresh_btn'
  ),
  '[-]Experimental' = list(
    'experimental',
    'limreal',
    'limimag'
  )
)

# End of input
# ----------------------------------  Outputs ----------------------------------
#' Define Outputs
#'
#' Use function `define_output` to define outputs.
#' Make sure rave toolbox (dev_[your_package_name]) is loaded.
#'
#' @Usage: define_output(definition, title, width, order)
#'
#' @param definition defines output types, for example:
#'     verbatimTextOutput('text_result')
#'   defines output type that dumps whatever is printed by function `text_result`
#'
#' Here are some possible types
#'   1. textOutput is an output for characters
#'   2. verbatimTextOutput is for console print
#'   3. plotOutput is for figures
#'   4. tableOutput is to tables
#'   5. customizedUI is for advanced UI controls
#'
#'   There are lots of output types and R packages such as DT, threejsr can provide
#'   very useful output types. Please check vignettes.
#'
#'   The easiest way to look for usage is using `help` function:
#'   help('verbatimTextOutput'), or ??verbatimTextOutput
#'
#'
#' @param title is the title for output
#'
#' @param width an integer from 1 to 12, defines the percentage of output width
#'   12 means 100% width, 6 means 50% and 4 means 33% width.
#'
#' @param order numeric order of outputs. Outputs will be re-ordered by this argument
#'

define_output(
  definition = verbatimTextOutput('current_sel'),
  title = 'Currently Loaded Trials',
  width = 5,
  order = 1
)

define_output(
  definition = verbatimTextOutput('possible_sel'),
  title = 'Available (Processed) Trials',
  width = 7,
  order = 2
)

define_output(
  plotOutput('fragility_map'),
  title = 'Fragility Map',
  width = 9,
  order = 3
)

define_output(
  definition = tableOutput('fragility_table'),
  title = 'Most/Least Fragile Electrodes',
  width = 3,
  order = 4
)

define_output(
  definition = plotOutput('voltage_trace', height = '80vh'),
  title = 'EEG Voltage Trace Viewer',
  width=12,
  order = 5
)

define_output_3d_viewer(
  outputId = 'power_3d',
  message = 'Click here to reload viewer',
  title = 'Results on surface',
  height = '500px',
  order = 6
)

# output_layout <- list(
#   'Status' = list(
#     'current_sel',
#     'possible_sel'
#   ),
#   'Visualization' = list(
#     'voltage_trace',
#     'fragility_map',
#     'fragility_table',
#     'power_3d'
#   )
# )

# <<<<<<<<<<<< End ----------------- [DO NOT EDIT THIS LINE] -------------------

# -------------------------------- View layout ---------------------------------

# Preview

view_layout('fragility')
