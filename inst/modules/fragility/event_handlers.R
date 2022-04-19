input = getDefaultReactiveInput()
output = getDefaultReactiveOutput()
session = getDefaultReactiveDomain()
static_data = dipsaus::fastmap2()

local_data = reactiveValues(
  v = NULL,
  check = NULL,
  pt_info = NULL,
  adj_info = NULL,
  brain_f = NULL,
  J = NULL,
  estimate = NULL,
  Hsize = NULL,
  selected = list(
    adj = '',
    f = ''
  )
)

observeEvent(
  input$process_pt, {
    print('process_pt button clicked')
    showNotification('Pre-processing patient...', id = 'pt_processing')
    local_data$pt_info <- process_fragility_patient(
      v = local_data$v,
      unit = input$recording_unit,
      srate = srate,
      halve = FALSE
    )
    saveRDS(local_data$pt_info, file = paste0(module_data,'/',subject_code,'_pt_info'))
    removeNotification(id = 'pt_processing')
    showNotification('Pre-processing complete!')
  }
)

observeEvent(
  input$gen_adj, {
    if (local_data$check$pt) {
      showNotification('Calculating estimated time...', id = 'loadingModal')
      # Show estimated adj array calculation time
      local_data$est <- estimate_time(local_data$pt_info, as.numeric(requested_tstep), requested_twindow)
      removeNotification(id = 'loadingModal')
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
          p('Required futures.globals.maxSize: '),
          tags$blockquote(format(local_data$est$Hsize, units = 'MB')),
          p('This might take a while.'),
          hr(),
          local_data$est$time
        ))
      ))
    } else {
      showNotification('Patient has not been processed yet. Please click the "Pre-process Patient" button under "Load Patient".', duration = 10)
    }
  }
)

observeEvent(input$cancel, {
  removeModal()
})

observeEvent(
  input$ok, {
    options(future.globals.maxSize = local_data$est$Hsize * 1024^2)
    local_data$adj_info <- generate_adj_array(
      t_window = input$requested_twindow,
      t_step = as.numeric(input$requested_tstep),
      v = local_data$pt_info$v,
      trial_num = tnum_adj,
      nlambda = input$requested_nlambda,
      ncores = input$requested_ncores
    )
    local_data$adj_info <- append(local_data$adj_info, list(trial = tnum_adj))
    saveRDS(local_data$adj_info, file = paste0(module_data,'/',subject_code,'_adj_info_trial_',tnum_adj))
    local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
    showNotification('Adjacency array generation finished!')
    removeModal()
  }
)

observeEvent(
  input$gen_f, {
    if (local_data$check$pt) {
      if (local_data$check$adj[tnum_adj]) {
        print('gen_f button clicked')
        f_info <- generate_fragility_matrix(
          A = adj_info$A,
          elec = attr(pt_info$v, "dimnames")$Electrode
        )
        f_info <- append(f_info, list(trial = tnum_adj))
        saveRDS(f_info, file = paste0(module_data,'/',subject_code,'_f_info_trial_',tnum_adj))
        local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
        updateSelectInput(session = session, inputId = 'requested_conditions',
                          choices = module_tools$get_meta('trials')$Condition[local_data$check$f],
                          selected = input$requested_conditions)
      } else {
        showNotification('This trial does not yet have an adjacency array, which is required for generating the fragility matrix. Click the "Generate Adjacency Array" button first.')
      }
    } else {
      showNotification('Patient has not been processed yet. Please click the "Pre-process Patient" button under "Load Patient".', duration = 10)
    }
  }
)

observeEvent(
  input$adj_conditions, {
    t <- trial$Trial[trial$Condition %in% input$adj_conditions]
    local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
    if (local_data$check$adj[t]) {
      if (is.null(local_data$adj_info)) {
        # if the file exists but hasn't been loaded in yet
        local_data$adj_info <- readRDS(paste0(module_data,'/',subject_code,'_adj_info_trial_',t))
      } else if (local_data$adj_info$trial != t) {
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
)

observeEvent(
  input$requested_conditions, {
    updateNumericInput(session = session, inputId = 'f_list_length', 
                       max = floor(length(preload_info$electrodes)/2), 
                       value = floor(length(preload_info$electrodes)/8))
  }
)

observeEvent(
  input$refresh_btn, {
    showNotification('Re-reading patient info...', id = 'refreshing')
    local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
    
    local_data$pt_info <- readRDS(paste0(module_data,'/',subject_code,'_pt_info'))
    
    updateSelectInput(session = session, inputId = 'requested_conditions',
                      choices = module_tools$get_meta('trials')$Condition[local_data$check$f],
                      selected = input$requested_conditions)
    removeNotification('refreshing')
  }
)

observeEvent(
  input$test, {
    print(str(local_data$pt_info))
  }
)

power_3d_fun = function(need_calc, side_width, daemon_env, viewer_proxy, ...){
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
    values = isolate(local_data$brain_f)
    
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

