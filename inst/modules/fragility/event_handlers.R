# event_handlers.R defines what code is run when different interactive events
# are triggered.

input = getDefaultReactiveInput()
output = getDefaultReactiveOutput()
session = getDefaultReactiveDomain()
static_data = dipsaus::fastmap2()

local_data = reactiveValues(
  v = NULL,
  check = NULL,
  pt_info = NULL,
  adj_info = NULL,
  f_info = NULL,
  requested_electrodes = NULL,
  v_electrodes = NULL,
  twindow = NULL,
  tstep = NULL,
  brain_f = NULL,
  est = NULL,
  J = NULL,
  v_loaded = NULL,
  vmat_params = NULL,
  f_plot_params = NULL,
  f_table_params = NULL,
  selected = list(
    adj = '',
    f = ''
  )
)

# Pre-process patient button clicked
observeEvent(
  input$process_pt, {
    print('process_pt button clicked')
    showNotification('Pre-processing patient...', id = 'pt_processing', duration = NULL)
    
    # volt <- module_tools$get_voltage()
    # local_data$v <- volt$get_data()
    
    # get processed pt_info
    local_data$pt_info <- process_fragility_patient(
      v = local_data$v,
      unit = input$recording_unit,
      srate = srate,
      halve = input$half_hz
    )
    
    # save into module_data folder
    saveRDS(local_data$pt_info, file = paste0(module_data,'/',subject_code,'_pt_info'))
    removeNotification(id = 'pt_processing')
    showNotification('Pre-processing complete!')
  }
)

# Generate adjacency array button clicked
observeEvent(
  input$gen_adj, {
    
    # re-check what files are available
    # local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
    
    if (local_data$check$pt) {
      showNotification('Estimating required futures.globals.maxSize...', id = 'loading_modal', duration = NULL)
      
      # convert requested_tstep and twindow from ms to # of datapoints within that timewindow using Hz
      local_data$twindow <- 2 * round(input$requested_twindow * local_data$pt_info$srate / 2000) # twindow needs to be even number
      if (as.numeric(input$requested_tstep) == input$requested_twindow) {
        local_data$tstep <- local_data$twindow
      } else if (as.numeric(requested_tstep)*2 == requested_twindow) {
        local_data$tstep <- local_data$twindow/2
      } else {
        stop('Time step should always be equal to or half of time window')
      }
      
      # calculate estimated time (currently working on this feature)
      local_data$est <- estimate_time(local_data$pt_info, local_data$tstep, local_data$twindow)
      
      removeNotification(id = 'loading_modal')
      
      # pop-up screen with additional info for adj array processing
      showModal(modalDialog(
        title = 'Confirmation',
        easyClose = F,
        footer = tagList(
          actionButton(inputId = ns('cancel'), 'Cancel'),
          actionButton(inputId = ns('ok'), 'Do it!')
        ),
        fluidRow(column(
          width = 12,
          p('Generate adjacency array for: '),
          tags$blockquote(paste0(local_data$est$J, ' time windows, ', length(preload_info$electrodes), ' electrodes (', dipsaus::deparse_svec(preload_info$electrodes), ')')),
          #p('Required futures.globals.maxSize: '),
          #tags$blockquote(format(local_data$est$Hsize, units = 'MB')),
          p('This might take a while. Estimated time remaining will appear after first time window is calculated.')
          # hr(),
          # local_data$est$time
        ))
      ))
    } else {
      showNotification('Patient has not been processed yet. Please click the "Pre-process Patient" button under "Load Patient".', duration = 10)
    }
  }
)

# cancel button on pop-up screen is clicked
observeEvent(input$cancel, {
  removeModal()
})

# "Do it!" button on pop-up screen is clicked (start adj array processing)
observeEvent(
  input$ok, {
    # set max future.globals export size to whatever is necessary for parallel processing
    options(future.globals.maxSize = local_data$est$Hsize * 1024^2)
    # options(future.globals.maxSize = 500 * 1024^2)
    
    # get adj_info
    local_data$adj_info <- generate_adj_array(
      t_window = local_data$twindow,
      t_step = local_data$tstep,
      v = local_data$pt_info$v,
      trial_num = tnum_adj,
      nlambda = input$requested_nlambda,
      ncores = input$requested_ncores
    )
    
    # add trial metadata and save adj_info with trial number into module_data
    local_data$adj_info <- append(local_data$adj_info, list(trial = tnum_adj))
    saveRDS(local_data$adj_info, file = paste0(module_data,'/',subject_code,'_adj_info_trial_',tnum_adj))
    showNotification('Adjacency array generation finished!')
    removeModal()
  }
)

# Generate fragility matrix button clicked
observeEvent(
  input$gen_f, {
    
    # re-check what files are available
    # local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
    
    # show different messages based on what files are available
    if (local_data$check$pt) {
      if (local_data$check$adj[tnum_adj]) {
        print('gen_f button clicked')
        local_data$f_info <- generate_fragility_matrix(
          A = local_data$adj_info$A,
          elec = attr(local_data$pt_info$v, "dimnames")$Electrode
        )
        local_data$f_info <- append(local_data$f_info, list(trial = tnum_adj))
        saveRDS(local_data$f_info, file = paste0(module_data,'/',subject_code,'_f_info_trial_',tnum_adj))
        # local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
        updateSelectInput(session = session, inputId = 'requested_conditions',
                          choices = module_tools$get_meta('trials')$Condition[local_data$check$f],
                          selected = input$requested_conditions)
      } else {
        showNotification('This trial does not yet have an adjacency array, which is required for generating the fragility matrix. Click the "Generate Adjacency Array" button first.')
      }
    } else {
      showNotification('Patient has not been processed yet. Please click the "Pre-process Patient" button under "Load Patient".')
    }
  }
)

observeEvent(
  input$adj_conditions, {
    if (exists('trial')) {
      t <- trial$Trial[trial$Condition %in% input$adj_conditions]
      # local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
      if (shiny::isTruthy(local_data$adj_info) & local_data$check$adj[t]) {
        if (local_data$adj_info$trial != t) {
          # if the user requests adj_info for a different trial
          local_data$adj_info <- readRDS(paste0(module_data,'/',subject_code,'_adj_info_trial_',t))
        }
        # save currently selected trial for "Currently Loaded Trials" display
        local_data$selected$adj <- local_data$adj_info$trial
      }
      updateActionButton(session = session, inputId = 'gen_adj', 
                         label = paste0('Generate Adjacency Array for ', input$adj_conditions))
      updateActionButton(session = session, inputId = 'gen_f', 
                         label = paste0('Generate Fragility Matrix for ', input$adj_conditions))
    }
  }
)

observeEvent(
  input$text_electrode, {
    local_data$requested_electrodes <- dipsaus::parse_svec(input$text_electrode)
    if (shiny::isTruthy(input$auto_calc) & !is.null(local_data$requested_electrodes)) {
      if (!all(local_data$requested_electrodes %in% preload_info$electrodes)) {
        stop('Please only select loaded electrodes.')
      }
      
      if (input$v_sync) {
        updateTextInput(session = session, inputId = 'v_electrode', value = input$text_electrode)
      }
      
      updateNumericInput(session = session, inputId = 'f_list_length', 
                         max = floor(length(local_data$requested_electrodes)/2), 
                         value = max(c(min(c(floor(length(local_data$requested_electrodes)/2), 10)), 1)))
      
      if (!is.null(input$requested_conditions)) {
        showNotification('Updating fragility map...', id = 'updating_f', duration = NULL)
        tnum <- trial$Trial[trial$Condition %in% input$requested_conditions]
        f_outputs <- draw_f_map_table(
          tnum = tnum, 
          adj_info = local_data$adj_info,
          f_path = paste0(module_data,'/',subject_code,'_f_info_trial_'), 
          subject_code = subject_code,
          requested_electrodes = local_data$requested_electrodes,
          sort_fmap = input$sort_fmap,
          check = local_data$check,
          f_list_length = input$f_list_length
        )
        
        if (input$exp_fmap){
          f_outputs$f_plot_params$mat <- f_outputs$f_plot_params$mat^input$exponentiate
        }
        
        local_data$brain_f <- f_outputs$brain_f
        local_data$f_plot_params <- f_outputs$f_plot_params
        local_data$f_table_params <- f_outputs$f_table_params
        local_data$selected$f <- f_outputs$sel
        
        local_data$J <- length(f_outputs$f_plot_params$x)
        local_data$f_plot_params <- append(local_data$f_plot_params, list(sz_onset = input$sz_onset))
        
        local_data$brain_f_plot <- dplyr::mutate(local_data$brain_f, Avg_Fragility = Avg_Fragility^input$exponentiate)
        rebuild_3d_viewer()
        
        removeNotification('updating_f')
      } else {
        local_data$brain_f <- NULL
        local_data$f_plot_params <- NULL
        local_data$f_table_params <- NULL
        local_data$selected$f <- ''
      }
    }
  }
)

observeEvent(
  list(
    input$requested_conditions,
    input$sort_fmap,
    input$exp_fmap,
    input$f_list_length,
    input$sz_onset
  ), {
    local_data$requested_electrodes <- dipsaus::parse_svec(input$text_electrode)
    if (shiny::isTruthy(input$auto_calc) & !is.null(local_data$requested_electrodes)) {
      if (!is.null(input$requested_conditions)) {
        showNotification('Updating fragility map...', id = 'updating_f', duration = NULL)
        tnum <- trial$Trial[trial$Condition %in% input$requested_conditions]
        f_outputs <- draw_f_map_table(
          tnum = tnum, 
          adj_info = local_data$adj_info,
          f_path = paste0(module_data,'/',subject_code,'_f_info_trial_'), 
          subject_code = subject_code,
          requested_electrodes = local_data$requested_electrodes,
          sort_fmap = input$sort_fmap,
          check = local_data$check,
          f_list_length = input$f_list_length
        )
        
        if (input$exp_fmap){
          f_outputs$f_plot_params$mat <- f_outputs$f_plot_params$mat^input$exponentiate
        }
        
        local_data$brain_f <- f_outputs$brain_f
        local_data$f_plot_params <- f_outputs$f_plot_params
        local_data$f_table_params <- f_outputs$f_table_params
        local_data$selected$f <- f_outputs$sel
        
        local_data$J <- length(f_outputs$f_plot_params$x)
        local_data$f_plot_params <- append(local_data$f_plot_params, list(sz_onset = input$sz_onset))
        
        local_data$brain_f_plot <- dplyr::mutate(local_data$brain_f, Avg_Fragility = Avg_Fragility^input$exponentiate)
        rebuild_3d_viewer()
        
        removeNotification('updating_f')
      } else {
        local_data$brain_f <- NULL
        local_data$f_plot_params <- NULL
        local_data$f_table_params <- NULL
        local_data$selected$f <- ''
      }
    }
  }
)

observeEvent(
  input$draw_f_map, {
    local_data$requested_electrodes <- dipsaus::parse_svec(input$text_electrode)
    if (!all(local_data$requested_electrodes %in% preload_info$electrodes)) {
      stop('Please only select loaded electrodes.')
    }
    showNotification('Updating fragility map...', id = 'updating_f', duration = NULL)
    tnum <- trial$Trial[trial$Condition %in% input$requested_conditions]
    f_outputs <- draw_f_map_table(
      tnum = tnum, 
      adj_info = local_data$adj_info,
      f_path = paste0(module_data,'/',subject_code,'_f_info_trial_'), 
      subject_code = subject_code,
      requested_electrodes = local_data$requested_electrodes,
      sort_fmap = input$sort_fmap,
      check = local_data$check,
      f_list_length = input$f_list_length
    )
    
    if (input$exp_fmap){
      f_outputs$f_plot_params$mat <- f_outputs$f_plot_params$mat^input$exponentiate
    }
    
    local_data$brain_f <- f_outputs$brain_f
    local_data$f_plot_params <- f_outputs$f_plot_params
    local_data$f_table_params <- f_outputs$f_table_params
    local_data$selected$f <- f_outputs$sel
    
    local_data$J <- length(f_outputs$f_plot_params$x)
    local_data$f_plot_params <- append(local_data$f_plot_params, list(sz_onset = input$sz_onset))
    
    local_data$brain_f_plot <- dplyr::mutate(local_data$brain_f, Avg_Fragility = Avg_Fragility^input$exponentiate)
    rebuild_3d_viewer()
    
    removeNotification('updating_f')
  }
)

observeEvent(
  input$exponentiate, {
    if (shiny::isTruthy(local_data$brain_f)) {
      local_data$brain_f_plot <- dplyr::mutate(local_data$brain_f, Avg_Fragility = Avg_Fragility^input$exponentiate)
      rebuild_3d_viewer()
    }
  }
)

observeEvent(
  input$load_v_traces, {
    t <- trial$Trial[trial$Condition %in% input$v_conditions]
    # local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
    local_data$voltage_electrodes <- dipsaus::parse_svec(input$v_electrode)
    if (!all(local_data$voltage_electrodes %in% preload_info$electrodes)) {
      stop('Please only select loaded electrodes.')
    }
    if (local_data$check$elist) {
      elec_labels <- paste0(local_data$check$elec_list$Label[local_data$voltage_electrodes], '(', local_data$voltage_electrodes, ')')
    } else {
      elec_labels <- local_data$voltage_electrodes
    }
    local_data$vmat_params <- list(
      mat = local_data$v[t,,as.character(local_data$voltage_electrodes)],
      elec_labels = elec_labels
    )
    local_data$v_loaded <- TRUE
  }
)

observeEvent(
  list(
    input$v_conditions,
    input$v_electrode,
    input$v_sync
  ), {
    if (exists('trial') & shiny::isTruthy(local_data$v_loaded)) {
      t <- trial$Trial[trial$Condition %in% input$v_conditions]
      # local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
      if (input$v_sync){
        updateTextInput(session = session, inputId = 'v_electrode', value = input$text_electrode)
        local_data$voltage_electrodes <- dipsaus::parse_svec(input$v_electrode)
      } else {
        local_data$voltage_electrodes <- dipsaus::parse_svec(input$v_electrode)
      }
      if (!all(local_data$voltage_electrodes %in% preload_info$electrodes)) {
        stop('Please only select loaded electrodes.')
      }
      if (local_data$check$elist) {
        elec_labels <- paste0(local_data$check$elec_list$Label[local_data$voltage_electrodes], '(', local_data$voltage_electrodes, ')')
      } else {
        elec_labels <- local_data$voltage_electrodes
      }
      local_data$vmat_params <- list(
        mat = local_data$v[t,,as.character(local_data$voltage_electrodes)],
        elec_labels = elec_labels
      )
    }
  }
)

observeEvent(
  input$refresh_btn, {
    showNotification('Re-reading patient info...', id = 'refreshing', duration = NULL)
    local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
    if (local_data$check$pt) {
      local_data$pt_info <- readRDS(paste0(module_data,'/',subject_code,'_pt_info'))
    }
    
    if (local_data$check$adj[tnum_adj]) {
      local_data$adj_info <- readRDS(paste0(module_data,'/',subject_code,'_adj_info_trial_',tnum_adj))
    }
    
    updateSelectInput(session = session, inputId = 'requested_conditions',
                      choices = module_tools$get_meta('trials')$Condition[local_data$check$f],
                      selected = input$requested_conditions)
    removeNotification('refreshing')
  }
)

# observeEvent(
#   input$test, {
#     local_data$brain_f$Avg_Fragility <- local_data$brain_f$Avg_Fragility^input$exponentiate
#     rebuild_3d_viewer()
#   }
# )

power_3d_fun <- function(need_calc, side_width, daemon_env, viewer_proxy, ...){
  showNotification(p('Rebuild 3d viewer...'), id='power_3d_fun')
  brain <- rave::rave_brain2(subject = subject)
  
  if(is.null(brain)){
    rave::close_tab('power_explorer', 'Results on surface')
    showNotification('No surface file is available...', duration=2)
  }
  
  shiny::validate(shiny::need(!is.null(brain), message = 'No surface/volume file found!'))
  re = NULL
  
  display_side = isTRUE(isolate(input$power_3d_widget_side_display))
  
  # zoom_level = shiny::isolate(viewer_proxy$main_camera$zoom)
  # controllers = viewer_proxy$get_controllers()
  # 
  # if(isolate(input$synch_3d_viewer_bg)) {
  #   bgcolor = isolate(input$background_plot_color_hint)
  #   if(bgcolor == 'gray') {
  #     bgcolor = '#1E1E1E'
  #   } else {
  #     bgcolor %<>% col2hex
  #   }
  #   controllers[['Background Color']] = bgcolor
  # }
  
  if( need_calc ){
    values = local_data$brain_f_plot
    
    # values = localdata.frame("Subject"=subject$subject_code,"Electrode"=requested_electrodes,"Time"=0,"Avg_Fragility"=local_data$avg_fragility)
    # values$Subject = as.factor(subject$subject_code)
    # values$Electrode = as.numeric(colnames(or))
    # values$Time = 0
    
    # if(!is.null(values$Passes_Filters)) {
    #   values$Passes_Filters[values$Passes_Filters==1] = 'Pass'
    #   values$Passes_Filters %<>% factor(levels=c('Pass'))
    # }
    # if(!is.null(values$Selected_Electrodes)) {
    #   values$Selected_Electrodes[values$Selected_Electrodes==1] = 'Selected'
    #   values$Selected_Electrodes %<>% factor(levels='Selected')
    # }
    
    # let's also set the electrode freesurfer label into the dset if we have it
    # if('freesurferlabel' %in% names(electrodes_csv)) {
    #   # I'm randomizing the factor order here so the colors will not be nearby
    #   values$Anatomy = factor(electrodes_csv[['freesurferlabel']], levels=sample(unique(electrodes_csv[['freesurferlabel']])))
    # }
    
    brain$set_electrode_values(values)
    # assign('omnibus_results', omnibus_results, globalenv())
    
    # check to see if we've udpated the dependent variable. We do this by comparing this list of Possible DVs with the 
    # actual current DV
    # curr_dv =  controllers[['Display Data']]
    # new_dv = rownames(or)[1]
    # is_old_dv = curr_dv %in% str_subset(format_unit_of_analysis_name(get_unit_of_analysis(names=T)),
    #                                     new_dv, negate = TRUE)
    # 
    # if(!length(controllers[['Display Data']]) || controllers[['Display Data']] == '[None]' || is_old_dv) {
    #   controllers[['Display Data']] = new_dv
    #   v = ceiling(max(abs(or[1,])) )
    #   controllers[['Display Range']] = sprintf('-%s,%s', v, v )
    #   # tr = controllers[['Threshold Range']]
    # }
    # 
    ### maybe we don't always want to show legend...
    # controllers[['Show Legend']] = TRUE
    
    ##FIXME -- change this to use the 3dViewer heatmap selector. not yet built
    # cp = input$viewer_color_palette
    # if(is.null(cp) || cp == 'Synch with heatmaps') {
    #   .colors = get_heatmap_palette(input$heatmap_color_palette)
    # } else {
    #   .colors = get_heatmap_palette(cp)
    # }
    # pal = expand_heatmap(.colors, ncolors=128)
    # pval_pal = expand_heatmap(
    #   rev(tail(.colors, ceiling(length(.colors)/2))),
    #   ncolors=128, bias=10)
    # pals = list(pal)
    # pals[2:nrow(omnibus_results )] = pals
    # names(pals) = fix_name_for_js(rownames(omnibus_results))
    # 
    # pals[str_detect(names(pals), 'p\\.')] = list(pval_pal)
    # pals[names(pals) %in% c('p')] = list(pval_pal)
    # 
    # pals[['[Subject]']] = 'black'
    # pals[['Passes_Filters']] = 'black'
    
    # re = brain$plot(symmetric = 0, # palettes = pals,
    #                 side_width = side_width / 2, side_canvas = TRUE, 
    #                 side_display = display_side, start_zoom = zoom_level, controllers = controllers,
    #                 control_presets = 'syncviewer', timestamp = FALSE)
    re = brain$plot()
    
  }else{
    # optional, if you want to change the way 3D viewer looks in additional tab
    daemon_env$widget = brain$plot(side_width = side_width, side_canvas = FALSE, side_display = display_side)
    
    # Just initialization, no need to show sidebar
    re = brain$plot(side_width = side_width / 2, side_canvas = TRUE, side_display = display_side,
                    control_presets = 'syncviewer', timestamp = FALSE, control_display=FALSE)
  }
  
  
  
  shiny::removeNotification('power_3d_fun')
  
  re
}

rebuild_3d_viewer <- function() {
  if(rave::rave_context()$context %in% c('rave_running', 'default')) {
    dipsaus::cat2("Updating 3D viewer from r3dv")
    btn_val = (input[['power_3d_btn']]) - 0.001
    dipsaus::set_shiny_input(session = session, inputId = 'power_3d_btn',
                             value = btn_val, method = 'proxy', priority = 'event')
  }
}

