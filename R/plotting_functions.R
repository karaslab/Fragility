#' draw_heat_map
#' this function makes an image using matrix mat
#' @param mat a 2D matrix
#' @param x vector for x axis
#' @param y vector for y axis
#' @param zlim vector of length 2 specifying range of image
#'
#' @return heatmap image
#' @export
#'
#' @examples
#' draw_heat_map(mat, x, y)
draw_heat_map <- function(mat, x, y, zlim) {
  image(z = t(mat), x = x, y = y,
        zlim = zlim,axes=F,ylab='frequency',xlab='time',
        col = colorRampPalette(c('blue', 'green', 'red'))(101))
}