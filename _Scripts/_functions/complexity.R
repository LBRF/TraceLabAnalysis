### Figure complexity analysis functions for TraceLab ###


### Bezier functions ###

line_length <- function(start.x, start.y, end.x, end.y) {
  dx <- (end.x - start.x)
  dy <- (end.y - start.y)
  sqrt(dx ** 2 + dy ** 2)
}


# Based on equation & C code from the following link:
# http://segfaultlabs.com/docs/quadratic-bezier-curve-length
#	(link dead, but available via wayback)
bezier_length <- function(start.x, start.y, end.x, end.y, ctrl.x, ctrl.y) {

  ax <- start.x - 2 * ctrl.x + end.x
  ay <- start.y - 2 * ctrl.y + end.y
  bx <- 2 * ctrl.x - 2 * start.x
  by <- 2 * ctrl.y - 2 * start.y

  A <- 4 * (ax ** 2 + ay ** 2)
  B <- 4 * (ax * bx + ay * by)
  C <- bx ** 2 + by ** 2

  sqrtABC <- 2 * sqrt(A + B + C)
  A2 <- sqrt(A)
  A32 <- 2 * A * A2
  C2 <- 2 * sqrt(C)
  BA <- B / A2

  suppressWarnings({
    n1 <- A32 * sqrtABC + A2 * B * (sqrtABC - C2)
    n2 <- ((4 * C * A) - B ** 2) * log((2 * A2 + BA + sqrtABC) / (BA + C2))
    len <- (n1 + n2) / (4 * A32)
  })

  # If length is NA or NaN, recalculate length of segment as straight line
  len <- ifelse(is.na(len), line_length(start.x, start.y, end.x, end.y), len)

  len
}


get_fig_points <- function(t, start.x, start.y, end.x, end.y, ctrl.x, ctrl.y) {

  # Define constants and transition values
  ax <- start.x - ctrl.x
  ay <- start.y - ctrl.y
  bx <- end.x - ctrl.x
  by <- end.y - ctrl.y

  # Create matrix for each segment's transitions
  tmat <- matrix(t, nrow = length(start.x), ncol = length(t), byrow = TRUE)

  pts <- tibble(
    x = as.vector(t(ctrl.x + ax * (1 - tmat) ** 2 + bx * tmat ** 2)),
    y = as.vector(t(ctrl.y + ay * (1 - tmat) ** 2 + by * tmat ** 2))
  )
  pts
}


bcurv <- function(t, start.x, start.y, end.x, end.y, ctrl.x, ctrl.y) {

  # Calculate first derivatives of bezier for x and y
  b1x <- (2 * (1 - t) * (ctrl.x - start.x)) + (2 * t * (end.x - ctrl.x))
  b1y <- (2 * (1 - t) * (ctrl.y - start.y)) + (2 * t * (end.y - ctrl.y))

  # Calculate second derivatives of bezier for x and y
  b2x <- (2 * (end.x - (2 * ctrl.x) + start.x))
  b2y <- (2 * (end.y - (2 * ctrl.y) + start.y))

  # Calculate signed curvature of bezier
  curv <- ((b1x * b2y) - (b1y * b2x)) / (((b1x^2) + (b1y^2)) ^ (3 / 2))

  curv
}



### Entropy functions ###

get_angle_diffs <- function(dx, dy, skip = 0) {

  if (skip > 0) {
    dx2 <- dx + lead(dx, skip)
    dy2 <- dy + lead(dy, skip)
  } else {
    dx2 <- dx
    dy2 <- dy
  }

  theta <- atan2(dy2, dx2) - atan2(lag(dy), lag(dx))
  theta <- ifelse(
    abs(theta) > pi,
    theta - (2 * pi) * sign(theta),
    theta
  )

  theta
}
