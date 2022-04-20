input = getDefaultReactiveInput()
output = getDefaultReactiveOutput()
session = getDefaultReactiveDomain()
static_data = dipsaus::fastmap2()

local_data = reactiveValues(
  subject_dir = module_tools$get_subject_dirs(),
  # module_data = module_tools$get_subject_dirs()$module_data_dir,
  # subject_code = subject$subject_code,
  # trial = module_tools$get_meta('trials'),
  # srate <- module_tools$get_sample_rate(original = TRUE),
  v = NULL,
  check = NULL,
  pt_info = NULL,
  adj_info = NULL,
  f_info = NULL,
  requested_electrodes = NULL,
  brain_f = NULL,
  J = NULL,
  estimate = NULL,
  Hsize = NULL,
  f_plot_params = NULL,
  f_table_params = NULL,
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
    local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
    if (local_data$check$pt) {
      showNotification('Calculating estimated time...', id = 'loading_modal')
      # Show estimated adj array calculation time
      local_data$est <- estimate_time(local_data$pt_info, as.numeric(requested_tstep), requested_twindow)
      removeNotification(id = 'loading_modal')
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
    showNotification('Adjacency array generation finished!')
    removeModal()
  }
)

observeEvent(
  input$gen_f, {
    local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
    if (local_data$check$pt) {
      if (local_data$check$adj[tnum_adj]) {
        print('gen_f button clicked')
        local_data$f_info <- generate_fragility_matrix(
          A = local_data$adj_info$A,
          elec = attr(local_data$pt_info$v, "dimnames")$Electrode
        )
        local_data$f_info <- append(local_data$f_info, list(trial = tnum_adj))
        saveRDS(local_data$f_info, file = paste0(module_data,'/',subject_code,'_f_info_trial_',tnum_adj))
        local_data$check <- check_subject(subject_code,subject_dir,trial$Trial)
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
      print('updating adj conditions')
      
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
  }
)

observeEvent(
  input$text_electrode, {
    local_data$requested_electrodes = dipsaus::parse_svec(text_electrode)
    updateNumericInput(session = session, inputId = 'f_list_length', 
                       max = floor(length(local_data$requested_electrodes)/2), 
                       value = max(c(min(c(floor(length(local_data$requested_electrodes)/2), 10)), 1)))
  }
)

observeEvent(
  list(
    input$requested_conditions,
    input$text_electrode,
    input$sort_fmap,
    input$f_list_length
  ), {
    if (shiny::isTruthy(input$auto_calc) & !is.null(local_data$requested_electrodes)) {
      local_data$requested_electrodes = dipsaus::parse_svec(text_electrode)
      if (!all(local_data$requested_electrodes %in% preload_info$electrodes)) {
        stop('Please only select loaded electrodes.')
      }
      if (!is.null(input$requested_conditions)) {
        showNotification('Updating fragility map...', id = 'updating_f')
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
        
        local_data$brain_f <- f_outputs$brain_f
        local_data$f_plot_params <- f_outputs$f_plot_params
        local_data$f_table_params <- f_outputs$f_table_params
        local_data$selected$f <- f_outputs$sel
        
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
# observeEvent(
#   input$requested_conditions, {
#     if (!is.null(local_data$requested_electrodes)) {
#       print('updating requested conditions')
#       showNotification('Updating fragility map...', id = 'updating_f')
#       tnum <- trial$Trial[trial$Condition %in% input$requested_conditions]
#       # print(tnum)
#       # if (length(tnum) > 1) {
#       #   f_plot <- list(
#       #     norm = matrix(data = 0, nrow = dim(local_data$adj_info$A)[1], ncol = dim(local_data$adj_info$A)[3]),
#       #     avg = vector(mode = 'numeric', length = dim(local_data$adj_info$A)[1]),
#       #     trial = numeric()
#       #   )
#       #   for (i in seq_along(tnum)) {
#       #     f_i <- readRDS(paste0(module_data,'/',subject_code,'_f_info_trial_',tnum[i]))
#       #     f_plot[1:2] <- mapply(function(x,y) x + y, f_plot[1:2], f_i[2:3])
#       #     f_plot$trial <- c(f_plot$trial,tnum[i])
#       #   }
#       #   f_plot[1:2] <- lapply(f_plot[1:2], function(x) x/length(tnum))
#       # } else {
#       #   f <- readRDS(paste0(module_data,'/',subject_code,'_f_info_trial_',tnum))
#       #   f_plot <- f[2:4]
#       # }
#       #
#       # f_plot$norm <- f_plot$norm[as.character(local_data$requested_electrodes),]
#       # f_plot$avg <- f_plot$avg[as.character(local_data$requested_electrodes)]
#       # local_data$brain_f <- data.frame("Subject"=subject_code,
#       #                                  "Electrode"=local_data$requested_electrodes,"Time"=0,
#       #                                  "Avg_Fragility"=f_plot$avg)
#       #
#       # elecsort <- sort(as.numeric(attr(f_plot$norm, "dimnames")[[1]]))
#       # fsort <- as.numeric(attr(sort(f_plot$avg), "names"))
#       #
#       # if (sort_fmap == 'Electrode') {
#       #   elec_order <- elecsort
#       # } else if (sort_fmap == 'Fragility') {
#       #   elec_order <- fsort
#       # }
#       #
#       # if (is.vector(f_plot$norm)){
#       #   elec_order <- local_data$requested_electrodes
#       #   x <- 1:length(f_plot$norm)
#       #   m <- t(t(f_plot$norm))
#       # } else {
#       #   f_plot$norm <- f_plot$norm[as.character(elec_order),]
#       #   x <- 1:dim(f_plot$norm)[2]
#       #   m <- t(f_plot$norm)
#       # }
#       #
#       # attr(m, 'xlab') = 'Time'
#       # attr(m, 'ylab') = 'Electrode'
#       # attr(m, 'zlab') = 'Fragility'
#       #
#       # if (local_data$check$elist) {
#       #   y <- paste0(local_data$check$elec_list$Label[elec_order], '(', elec_order, ')')
#       #   f_list <- paste0(local_data$check$elec_list$Label[fsort], '(', fsort, ')')
#       # } else {
#       #   y <- elec_order
#       #   f_list <- fsort
#       # }
#       #
#       # local_data$f_plot_params <- list(
#       #   mat = m,
#       #   x = x,
#       #   y = y,
#       #   zlim = c(0,1)
#       # )
#       #
#       # local_data$f_table_params <- data.frame(
#       #   # Ranking = 1:f_list_length,
#       #   Most.Fragile = rev(f_list)[1:f_list_length],
#       #   Least.Fragile = f_list[1:f_list_length]
#       # )
#       #
#       #
#       #
#       # local_data$selected$f <- f_plot$trial
#       f_outputs <- draw_f_map_table(
#         tnum = tnum,
#         adj_info = local_data$adj_info,
#         f_path = paste0(module_data,'/',subject_code,'_f_info_trial_'),
#         requested_electrodes = requested_electrodes
#       )
# 
#       local_data$brain_f <- f_outputs$brain_f
#       local_data$f_plot_params <- f_outputs$f_plot_params
#       local_data$f_table_params <- f_outputs$f_table_params
#       local_data$selected$f <- f_outputs$sel
# 
#       removeNotification('updating_f')
#     }
#   }
# )

observeEvent(
  input$draw_f_map, {
    showNotification('Updating fragility map...', id = 'updating_f')
    local_data$requested_electrodes = dipsaus::parse_svec(text_electrode)
    if (!all(local_data$requested_electrodes %in% preload_info$electrodes)) {
      stop('Please only select loaded electrodes.')
    }
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
    
    local_data$brain_f <- f_outputs$brain_f
    local_data$f_plot_params <- f_outputs$f_plot_params
    local_data$f_table_params <- f_outputs$f_table_params
    local_data$selected$f <- f_outputs$sel
    
    removeNotification('updating_f')
  }
)

observeEvent(
  input$refresh_btn, {
    showNotification('Re-reading patient info...', id = 'refreshing')
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

observeEvent(
  input$test, {
    print(module_data)
    print(subject_code)
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

