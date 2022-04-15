# Function names should be consistent with the output IDs
# will be used to render outputs

current_sel <- function(result, ...){
  sel <- result$get_value('selected')
  sel$f <- toString(sel$f)
  paste0('Adjacency Matrix: ', sel$adj, ' | Fragility Map: ', sel$f)
}

fragility_table <- function(result, ...) {
  f_table <- result$get_value('f_table_params')
  shiny::validate(shiny::need(!is.null(f_table), message = 'No valid fragility matrices detected!'))
  f_table
}

# least_fragile <- function(result, ...) {
#   f_text <- result$get_value('f_text_params')
#   if(is.null(f_text)) {
#     return('No valid fragility matrices detected! Please follow the steps on the left to generate one.')
#   } else {
#     paste0(f_text$num, ' least fragile electrodes: ', paste0(f_text$least, collapse = ', '))
#   }
# }

fragility_map <- function(result, ...) {
  f <- result$get_value('f_plot_params')
  shiny::validate(shiny::need(!is.null(f), message = 'No valid fragility matrices detected! Please follow the steps on the left to generate one.'))
  
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
    
    mtext(y, side=2, line=-1, at=yi, cex=(ravebuiltins:::rave_cex.lab*0.8), las=1)
  }, ravebuiltins:::spectrogram_heatmap_decorator())
  )
}
