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
  subject_dir <- module_tools$get_subject_dirs()$module_data_dir
  if (!file.exists(subject_dir)){
    dir.create(subject_dir)
  }
  subject_code <- subject$subject_code
  
  trial = module_tools$get_meta('trials')
  srate <- module_tools$get_sample_rate(original = TRUE)
  
  # check if subject has pt_info, adj_info, and f_info files
  check <- check_subject(subject_code,subject_dir,trial$Trial)
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
#'
define_input(
  definition = textInput(inputId = 'text_electrode', label = 'Electrode'),

  # Use help('textInput') to see the definition, or help('updateTextInput')
  # to see which input element(s) can be changed. In this case, they are:
  #   label, value, placeholder

  # We will change label and value
  # label will tell users which electrodes are loaded
  # value will be the first electrode

  init_args = c('label', 'value'),
  init_expr = {

    # check ?rave_prepare for what kind of data are loaded
    loaded_electrodes = preload_info$electrodes

    # Generate text for loaded electrodes
    text = dipsaus::deparse_svec(loaded_electrodes)

    # Update
    label = paste0('Select Electrode (Loaded: ', text, ')')
    value = loaded_electrodes
  }
)

define_input(
  definition = selectInput(inputId = 'recording_unit', label='Units of Recording', choices = character(0)),

  init_args = c('choices','selected'),

  init_expr = {
    choices = c('mV','uV','nV')
    selected = 'uV'
  }
)

define_input(
  selectInput(inputId = 'adj_conditions', choices = character(0), label='Seizure Trial for Matrix Generation'),

  init_args = c('choices', 'selected'),

  init_expr = {
    # the module_tools object allows access to the meta data, including the "trial label" variable called Condition
    # The trial numbers and condition labels, are in the epoch file
    # rstudioapi::navigateToFile(file.path(rave::rave_options('data_dir'), 'demo','DemoSubject','rave', 'meta','epoch_auditory_onset.csv'))

    choices = unique(module_tools$get_meta('trials')$Condition)
    selected = module_tools$get_meta('trials')$Condition[1]
  }
)

define_input(
  selectInput(inputId = 'requested_conditions', choices = character(0), multiple = TRUE, label='Seizure Trial(s) for Fragility Map Display'),

  init_args = c('label','choices', 'selected'),

  init_expr = {
    if (any(check$f)){
      choices = module_tools$get_meta('trials')$Condition[check$f]
      selected = module_tools$get_meta('trials')$Condition[check$f][1]
    } else {
      label = 'No valid fragility matrices detected! Please generate one above.'
    }

    # choices = module_tools$get_meta('trials')$Condition
    # selected = module_tools$get_meta('trials')$Condition[1]
  }
)

define_input(
  sliderInput(inputId = 'requested_twindow', label='Time Window Size', min=100, max=1000, value=250, step=50, ticks = FALSE),

  # could code in something to auto calculate timewindow using # of time points and sample rate?
  # init_args = c('min', 'max', 'value', 'step'),
  # 
  # init_expr = {
  #   # timepoints = length(preload_info$time_points)
  #   # hz = module_tools$get_sample_rate(original = TRUE)
  # }
)

define_input(
  selectInput(inputId = 'requested_tstep', label='Time Step', choices = character(0)),

  init_args = c('choices', 'selected'),

  init_expr = {
    choices = c(input$requested_twindow, input$requested_twindow/2)
    selected = input$requested_twindow/2
  }
)

define_input(
  numericInput(inputId = 'requested_nlambda', label='Number of Lambdas', min=1, max=100, value=16, step=1),

  # init_args = 'value',
  #
  # init_expr = {
  #   value =
  # }
)

define_input(
  numericInput(inputId = 'requested_ncores', label='Number of cores', min=1, max=32, value=8, step=1),

  # init_args = 'max',
  #
  # init_expr = {
  #   max = how to find max cores in RAVE settings?
  # }
)

define_input(
  numericInput(inputId = 'future_maxsize', label='future.globals.maxSize (in megabytes)', min=500, max=8000, value=2000, step=500),
)

define_input(
  selectInput(inputId = 'sort_fmap', choices = c('Electrode','Fragility'), selected = 'Electrode', label='Sort Fragility Map By:'),
)

define_input(
  actionButton(inputId = 'process_pt', label='Pre-process Patient'),
)

define_input(
  actionButton(inputId = 'gen_adj', label='Generate Adjacency Matrix for this Trial'),
)

define_input(
  actionButton(inputId = 'gen_f', label='Generate Fragility Matrix for this Trial'),
)

define_input(
  actionButton(inputId = 'refresh_btn', label='Refresh'),
)

# the input_layout list is used by rave to determine order and grouping of layouts
input_layout <- list(
  '[-]Load Patient' = list(
    'recording_unit',
    'process_pt'
  ),
  '[-]Adjacency Matrix' = list(
    'requested_twindow',
    'requested_tstep',
    'requested_nlambda',
    'requested_ncores',
    'future_maxsize',
    'adj_conditions',
    'gen_adj',
    'gen_f'
  ),
  'Fragility Map' = list(
    'requested_conditions',
    'text_electrode',
    'sort_fmap',
    'refresh_btn'
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
  width = 6,
  order = 1
)

define_output(
  plotOutput('fragility_map'),
  title = 'Fragility Map',
  width = 12,
  order = 2
)

#for some reason this makes main execute an additional time
define_output_3d_viewer(
  outputId = 'power_3d',
  message = 'Click here to reload viewer',
  title = 'Results on surface',
  height = '500px',
  order = 3,
)

# <<<<<<<<<<<< End ----------------- [DO NOT EDIT THIS LINE] -------------------

# -------------------------------- View layout ---------------------------------

# Preview

view_layout('fragility')
