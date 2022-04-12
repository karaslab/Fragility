input = getDefaultReactiveInput()
output = getDefaultReactiveOutput()
session = getDefaultReactiveDomain()
static_data = dipsaus::fastmap2()

local_data = reactiveValues(
  brain_f = NULL,
  J = NULL,
  estimate = NULL
)

power_3d_fun = function(need_calc, side_width, daemon_env, viewer_proxy, ...){
  showNotification(p('Rebuild 3d viewer...'), id='power_3d_fun')
  brain = rave::rave_brain2(subject = subject)
  
  if(is.null(brain)){
    rave::close_tab('power_explorer', 'Results on surface')
    showNotification('No surface file is available...', duration=2)
  }
  
  shiny::validate(shiny::need(!is.null(brain), message = 'No surface/volume file found!'))
  re = NULL
  
  display_side = isTRUE(isolate(input$power_3d_widget_side_display))
  # 
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
    values = local_data$brain_f
    
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

observeEvent(
  input$process_pt, {
    print('process_pt button clicked')
    volt <- module_tools$get_voltage()
    v <- volt$get_data()
    pt_info_all <- process_fragility_patient(
      v = v,
      unit = input$recording_unit,
      srate = srate,
      halve = FALSE
    )
    saveRDS(pt_info_all, file = paste0(subject_dir,'/',subject_code,'_pt_info'))
    print(trial)
    print(tnum)
    print(requested_electrodes)
  }
)

observeEvent(
  input$gen_adj, {
    showModal(
      modalDialog(
        title = 'Confirmation',
        easyClose = F,
        footer = tagList(
          actionButton(inputId = ns('cancel'), 'Cancel'),
          actionButton(inputId = ns('ok'), 'Do it!')
        ),
        fluidRow(
          column(
            width = 12,
            p('Generate adjacency array for: '),
            tags$blockquote(paste0(local_data$J, ' time windows, ', length(preload_info$electrodes), ' electrodes (', dipsaus::deparse_svec(preload_info$electrodes), ')')),
            p('This might take a while.'),
            hr(),
            local_data$estimate
          )
        )
      )
    )
  }
)

observeEvent(input$cancel, {
  # little trick: if server is running other tasks, modal will not be removed
  removeModal()
})

observeEvent(
  input$ok, {
    if (check$pt) {
      print('ok button clicked')
      print(input$requested_twindow)
      print(input$requested_tstep)
      print(str(pt_info_all$v))
      print(tnum_adj)
      print(input$requested_nlambda)
      print(input$requested_ncores)
      options(future.globals.maxSize = input$future_maxsize * 1024^2)
      adj_info <- generate_adj_array(
        t_window = input$requested_twindow,
        t_step = as.numeric(input$requested_tstep),
        v = pt_info_all$v,
        trial_num = tnum_adj,
        nlambda = input$requested_nlambda,
        ncores = input$requested_ncores
      )
      adj_info <- append(adj_info, list(trial = tnum_adj))
      saveRDS(adj_info, file = paste0(subject_dir,'/',subject_code,'_adj_info_trial_',tnum_adj))
    } else {
      showNotification('Patient has not been processed yet. Please click the "Pre-process Patient" button under "Load Patient".', duration = 10)
    }
  }
)

observeEvent(
  input$gen_f, {
    if (check$pt) {
      if (check$adj[tnum_adj]) {
        print('gen_f button clicked')
        f_info <- generate_fragility_matrix(
          A = adj_info$A,
          elec = attr(pt_info_all$v, "dimnames")$Electrode
        )
        f_info <- append(f_info, list(trial = tnum_adj))
        saveRDS(f_info, file = paste0(subject_dir,'/',subject_code,'_f_info_trial_',tnum_adj))
        check <- check_subject(subject$subject_code,module_tools$get_subject_dirs()$module_data_dir,trial$Trial)
        updateSelectInput(session = session, inputId = 'requested_conditions',
                          choices = module_tools$get_meta('trials')$Condition[check$f],
                          selected = input$requested_conditions)
      } else {
        showNotification('No valid adjacency arrays detected. Please choose a trial and click the "Generate Adjacency Matrix" button.', duration = 10)
      }
    } else {
      showNotification('Patient has not been processed yet. Please click the "Pre-process Patient" button under "Load Patient".', duration = 10)
    }
  }
)

observeEvent(
  input$adj_conditions, {
    print('update')
    updateActionButton(session = session, inputId = 'gen_adj', 
                       label = paste0('Generate Adjacency Array for ', input$adj_conditions))
    updateActionButton(session = session, inputId = 'gen_f', 
                       label = paste0('Generate Fragility Matrix for ', input$adj_conditions))
  }
)

observeEvent(
  input$refresh_btn, {
    check <- check_subject(subject$subject_code,module_tools$get_subject_dirs()$module_data_dir,trial$Trial)
    updateSelectInput(session = session, inputId = 'requested_conditions',
                      choices = module_tools$get_meta('trials')$Condition[check$f],
                      selected = input$requested_conditions)
  }
)
