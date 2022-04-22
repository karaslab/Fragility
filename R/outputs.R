# Function names should be consistent with the output IDs
# will be used to render outputs

current_sel <- function(result, ...){
  sel <- result$get_value('local_data')$selected
  if (sel$adj == ''){
    sel$adj <- 'None'
  }
  if (all(sel$f == '')){
    sel$f <- 'None'
  } else {
    sel$f <- toString(sel$f)
  }
  paste0('Adjacency Array: ', sel$adj, ' | Fragility Map: ', sel$f)
}

possible_sel <- function(result, ...){
  pos <- lapply(result$get_value('possible'), toString)
  if (all(pos$adj == '')){
    pos$adj <- 'None'
  }
  if (all(pos$f == '')){
    pos$f <- 'None'
  }
  paste0('Pt. Processed? ', pos$pt, ' | Adj. Array: ', pos$adj, ' | Frag. Map: ', pos$f)
}

fragility_table <- function(result, ...) {
  f_table <- result$get_value('local_data')$f_table_params
  shiny::validate(shiny::need(!is.null(f_table), message = 'No fragility map currently loaded!'))
  f_table
}

fragility_map <- function(result, ...) {
  f <- result$get_value('local_data')$f_plot_params
  shiny::validate(shiny::need(!is.null(f), message = 'No fragility map currently loaded! Please follow the steps on the left.'))
  
  y=f$y
  yi = seq_along(y)
  if(length(f$y) > 10) {
    .seq = seq(1, length(f$y), length.out=10)
    y = f$y[.seq]
    yi = .seq
  }
  
  ravebuiltins:::draw_many_heat_maps(list(
    list(
      data = f$mat,
      x = f$x,
      y = seq_along(f$y),
      has_trials = TRUE,
      range = 0:1
    )
  ), axes = c(TRUE,FALSE), PANEL.LAST = ravebuiltins:::add_decorator(function(...) {
    abline(v = f$sz_onset, lty = 2, lwd = 2)
    mtext(y, side=2, line=-1, at=yi, cex=(ravebuiltins:::rave_cex.lab*0.8), las=1)
  }, ravebuiltins:::spectrogram_heatmap_decorator())
  )
}
