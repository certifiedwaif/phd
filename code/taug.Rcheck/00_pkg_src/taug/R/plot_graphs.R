#' Plot graphs of the log-lower bound of tau_g for various values of n, p and R2
#'
#' @export
plot_graphs <- function()
{
  flag <- FALSE
  n <- 110
  R2 <- seq(0, 1, .01)
  cols <- colorRampPalette(c("blue", "red"), interpolate = "spline")(length(R2))
  for (p in 2:95) {
    y <- rep(NA, length(R2))
    for (i in seq_along(R2)) {
      y[i] <- taug::tau_g(n, p, R2[i])
    }
    if (!flag) {
      plot(R2, y, type = "l", col = cols[p], xlim = c(0, 1), ylim = c(-200, 100))
    } else {
        lines(R2, y, col = cols[p])
    }
    flag <- TRUE
  }
}

plot_matrix <- function()
{
  n <- 100
  R2 <- seq(0, .99, .01)
  cols <- colorRampPalette(c("blue", "red"), interpolate = "spline")(length(R2))
  p <- 2:95
  z <- matrix(NA, length(p), length(R2))
  for (p2 in p) {
    for (i in seq_along(R2)) {
      z[p2-1, i] <- taug::tau_g(n, p2, R2[i])
    }
  }
  library(Matrix)
  image(Matrix(z))
}

plot_surface <- function()
{
  n <- 100
  p <- 2:95
  R2 <- seq(0, .9, length.out=length(p))
  cols <- colorRampPalette(c("blue", "red"), interpolate = "spline")(length(R2))
  z <- matrix(NA, length(p), length(R2))
  for (p2 in p) {
    for (i in seq_along(R2)) {
      z[p2-1, i] <- taug::tau_g(n, p2, R2[i])
    }
  }
  library(rgl)
  rgl.surface(p, 100*R2, .5*z, color=cols, coords=c(1, 3, 2))
  title3d(main = "tau_g")
  axes3d()
  #rgl.close()
  #rgl.surface(1:100, 1:100, rep(100, 100^2))
}
