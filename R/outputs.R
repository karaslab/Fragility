# Function names should be consistent with the output IDs
# will be used to render outputs

estim_time <- function(result, ...){
  result$get_value('estimate')
}

current_sel <- function(result, ...){
  sel <- result$get_value('selected')
  sel$f <- toString(sel$f)
  paste0('Adjacency Matrix: ', sel$adj, '; Fragility Map: ', sel$f)
}

fragility_map <- function(result, ...) {
  f <- result$get_value('f_plot_params')
  
  draw_heat_map(mat = f$mat, x = f$x, y = f$y, zlim = f$zlim)
}
