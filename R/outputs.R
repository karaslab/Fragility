# Function names should be consistent with the output IDs
# will be used to render outputs

text_result <- function(result, ...){
  result$get_value('my_text_result')
}


fragility_map <- function(result, ...) {
  f <- result$get_value('f_plot_params')
  
  draw_heat_map(mat = f$mat, x = f$x, y = f$y, zlim = f$zlim)
}
