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
  cat(paste0('Adjacency Array: ', sel$adj), paste0('Fragility Map: ', sel$f), sep = '\n')
}

possible_sel <- function(result, ...){
  pos <- lapply(result$get_value('possible'), toString)
  if (all(pos$adj == '')){
    pos$adj <- 'None'
  }
  if (all(pos$f == '')){
    pos$f <- 'None'
  }
  cat(paste0('Pt. Processed? ', pos$pt), paste0('Adj. Array: ', pos$adj), paste0('Frag. Map: ', pos$f), sep = '\n')
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
  
  # map x axis from timewindows (f$x) to time (for mtext)
  xtime <- round(seq(f$tp[1], f$tp[length(f$tp)], length.out = 9), digits = 2)
  xi <- seq(1, length(f$x), length.out = 9)
  
  # map seizure onset from time (from slider input) to timewindows (for abline)
  secs <- seq(f$tp[1], f$tp[length(f$tp)])
  onset <- seq(1, length(f$x), length.out = length(secs))[match(f$sz_onset,secs)]
  
  ravebuiltins:::draw_many_heat_maps(list(
    list(
      data = f$mat,
      x = f$x,
      y = seq_along(f$y),
      has_trials = TRUE,
      range = 0:1
    )
  ), axes = c(FALSE,FALSE), PANEL.LAST = ravebuiltins:::add_decorator(function(...) {
    abline(v = onset, lty = 2, lwd = 2)
    mtext(y, side=2, line=-1, at=yi, cex=(ravebuiltins:::rave_cex.lab*0.8), las=1)
    mtext(xtime, side=1, line=1, at=xi, cex=(ravebuiltins:::rave_cex.lab*0.8), las=1)
  }, ravebuiltins:::spectrogram_heatmap_decorator())
  )
}

voltage_trace <- function(result, ...) {
  vmat_params <- result$get_value('local_data')$vmat_params
  shiny::validate(shiny::need(!is.null(vmat_params), message = 'Load the original EEG data with the button on the left.'))
  
  shiny::showNotification("Plotting voltage traces...", id = 'plot')
  
  # compress by factor of 4 to save plotting time
  vmat <- vmat_params$mat[seq(1,nrow(vmat_params$mat),4),]

  t <- dim(vmat)[1]*4
  N <- dim(vmat)[2]
  z <- (vmat - mean(vmat)) / sd(vmat)

  par(mar=c(3,6,0,0))
  rutabaga::plot_clean(as.numeric(attr(vmat, "dimnames")$Time), -3:N*3)
  axis(1)
  abline(h=0:(N-1)*3, col = "skyblue")
  for(ii in 1:ncol(vmat)) {
    lines(x = attr(vmat, "dimnames")$Time, y = z[,ii]+(ii-1)*3)
  }

  axis(2, at=0:(N-1)*3, vmat_params$elec_labels, las=1, tcl=0)
  
  shiny::removeNotification('plot')
}
