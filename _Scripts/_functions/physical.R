#### Physical tracing preprocessing and analysis functions for TraceLab ####


#### Import required packages ####

library(tibble)
library(vegan)
library(dtw)



#### Utility functions ####

# Calculates standard deviation of paired differences

paired.sd <- function(x) {
  SS <- sum(x ** 2)
  df <- length(x) - 1
  sqrt(SS / df)
}

# Reinterpolates a figure or tracing, either by time or by distance

reinterpolate <- function(x, y, time, equidistant = FALSE, fps = 60, n = NA) {

  # If equidistant, interpolate constant-velocity points along path. Otherwise,
  # use timestamps to generate points at same variable speed as input tracing.
  if (equidistant) {
    dists <- sqrt((x - lag(x))^2 + (y - lag(y))^2)
    dists[1] <- 0
    steps <- cumsum(dists)
  } else {
    steps <- time
  }

  # Figure out number of frames to generate based on tracing time and framerate,
  # unless number of frames is explicitly provided
  if (!is.na(n)) {
    out_num <- n
  } else {
    out_num <- round(max(time) / (1 / fps)) + 1 # add 1 because first time is 0
  }

  # Return new interpolated points
  interpolated <- tibble(
    x = approx(steps, x, n = out_num, method = "linear")$y,
    y = approx(steps, y, n = out_num, method = "linear")$y
  )

  interpolated
}

# Matches length of stimulus to response via reinterpolation

match_lengths <- function(x, y, tx, ty, equidist = FALSE) {

  # Get path lengths of stim and tracing
  stim_n <- max(which(!is.na(x)))
  trace_n <- max(which(!is.na(tx)))

  # Reinterpolate stim to make it equal to tracing in length and calculate
  # turn angles
  stime <- seq(1, stim_n) * (1 / 60)
  newstim <- reinterpolate(x[!is.na(x)], y[!is.na(y)], stime, n = trace_n)
  matched_df <- tibble(
    x = newstim$x, y = newstim$y,
    trace.x = tx[!is.na(tx)], trace.y = ty[!is.na(ty)],
    delta = get_angle_diffs(newstim$x, newstim$y),
    trace.delta = get_angle_diffs(tx[!is.na(tx)], ty[!is.na(ty)])
  )

  # If enabled, reinterpolate tracing so all points are equidistant and update
  # turn angles
  if (equidist) {
    ttime <- seq(1, trace_n) * (1 / 60)
    newtrace <- reinterpolate(
      tx[!is.na(tx)], ty[!is.na(ty)],
      ttime, equidistant = TRUE, n = trace_n
    )
    matched_df$trace.x <- newtrace$x
    matched_df$trace.y <- newtrace$y
    matched_df$trace.delta = get_angle_diffs(newtrace$x, newtrace$y)
  }

  matched_df
}



#### Transformation functions ####

procrustes2df <- function(x, y, tx, ty, delta) {

  stim <- cbind(x, y)
  resp <- cbind(tx, ty)

  # Do procrustes transformation and get translated dataframe
  proc <- procrustes(stim, resp, scale = TRUE, symmetric = FALSE)
  proc_mat <- predict(proc, resp)
  proc_df <- as_tibble(proc_mat, .name_repair = ~ c("trace.x", "trace.y"))

  # Bind stim columns to procrustes data
  proc_df <- add_column(proc_df, x = x, .before = 1)
  proc_df <- add_column(proc_df, y = y, .before = 2)
  proc_df <- add_column(proc_df, delta = delta, .before = 3)
  proc_df <- add_column(proc_df, trace.delta = get_angle_diffs(proc_df$trace.x, 
                                                  proc_df$trace.y))

  # Get procrustes metrics and add to dataframe
  proc_dx <- proc$xmean[1] - mean(tx)
  proc_dy <- proc$xmean[2] - mean(ty)
  proc_df$translation <- sqrt(proc_dx ** 2 + proc_dy ** 2)
  proc_df$scale <- proc$scale
  proc_df$rotation <- acos(proc$rotation[1, 1])

  proc_df
}

dtw2df <- function(x, y, tx, ty, step = mori2006, open_ends = FALSE, 
                   turn_angle=FALSE, absolute = FALSE, upsample = 1) {
  
  if (upsample > 1) {
    tmp <- reinterpolate(x, y, NA, TRUE, n = round(length(x) * upsample))
    x <- tmp$x
    y <- tmp$y
  }
  
  time <- 1:length(x) * (1 / 60)
  delta <- get_angle_diffs(x, y)
  trace.delta <- get_angle_diffs(tx, ty)
  
  if (turn_angle){
    d1 <- scale(trace.delta[!is.na(trace.delta)])
    d2 <- scale(delta[!is.na(delta)])
    if (absolute) {
      d1 <- abs(d1)
      d2 <- abs(d2)
    }
  } else {
    d1 <- cbind(tx[!is.na(tx)], ty[!is.na(ty)])
    d2 <- cbind(x[!is.na(x)], y[!is.na(y)])
  }
  
  tryCatch({
    warped <- dtw(
      x = d1,
      y = d2,
      step.pattern = step,
      window.type = "none",
      keep.internals = TRUE,
      distance.only = FALSE,
      open.end = open_ends,
      open.begin = open_ends
    )
    
    # Return dataframe with warped stimulus and response points
    dtw_df <- tibble(
      t_w = time[warped$index2],
      x = x[warped$index2],
      y = y[warped$index2],
      delta = delta[warped$index2],
      trace.t_w = time[warped$index1],
      trace.x = tx[warped$index1],
      trace.y = ty[warped$index1],
      trace.delta = trace.delta[warped$index1]
    ) %>%
      mutate(time = 1:n() * (1 / 60))
    
    dtw_df
  },
  error = function(x) {
    tibble(
      t_w = numeric(0),
      x = numeric(0),
      y = numeric(0),
      delta = numeric(0),
      trace.t_w = numeric(0),
      trace.x = numeric(0),
      trace.y = numeric(0),
      trace.delta = numeric(0)
    )
  })
}
