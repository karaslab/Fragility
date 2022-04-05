input = getDefaultReactiveInput()
output = getDefaultReactiveOutput()
session = getDefaultReactiveDomain()
static_data = dipsaus::fastmap2()


local_data = reactiveValues(
  
)

observeEvent(
  input$load_pt, {
    print('load_pt button clicked')
    volt <- module_tools$get_voltage()
    v <- volt$get_data()
    pt_info_all <- load_fragility_patient(
      v = v,
      unit = input$recording_unit,
      halve = FALSE
    )
    saveRDS(pt_info_all, file = paste0(subject_dir,'/',subject_code,'_pt_info'))
    print(trial)
    print(tnum)
    print(requested_electrodes)
  }
)

observeEvent(
  input$load_adj, {
    print('load_adj button clicked')
    print(input$requested_twindow)
    print(input$requested_tstep)
    print(str(pt_info_all$v))
    print(tnum_adj)
    print(input$requested_nlambda)
    print(input$requested_ncores)
    options(future.globals.maxSize = input$future_maxsize * 1024^2)
    adj_info <- generate_adj_array(
      t_window = input$requested_twindow,
      t_step = input$requested_tstep,
      v = pt_info_all$v,
      trial_num = tnum_adj,
      nlambda = input$requested_nlambda,
      ncores = input$requested_ncores
    )
    saveRDS(adj_info, file = paste0(subject_dir,'/',subject_code,'_adj_info_trial_',tnum_adj))
  }
)

observeEvent(
  input$load_f, {
    print('load_f button clicked')
    f_info <- generate_fragility_matrix(
      A = adj_info$A,
      elec = requested_electrodes
    )
    saveRDS(f_info, file = paste0(subject_dir,'/',subject_code,'_f_info_trial_',tnum))
  }
)


observeEvent(
  input$test, {
    print(input$recording_unit)
    print(input$requested_twindow)
  }
)
