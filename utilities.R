#   Copyright (C) 2017 Pierre-André MAUGIS
#   This code comes with ABSOLUTELY NO WARRANTY.
#   This is free software, and you are welcome to redistribute it
#   under certain conditions; see attached LICENSE file for details.
#
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or (at
#     your option) any later version.
#
#     This program is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
#     02110-1301 USA.
#
#
#

library(optimx)

ffct <- function(x, k) {
  # Falling factorial
  sapply(x, function(y) {
    prod(y:(y - k + 1))
  })
}

Schoose2 <- function(M) {
  # choose(M,2) for sparse matrix
  return(M * (M - 1) / 2)
}

Schoose3 <- function(M) {
  # choose(M,3) for sparse matrix
  return(M * (M - 1) * (M - 2) / 6)
}

emb_contour <- function(X, prec = 1e3, spar = .55) {
  ######
  # Input:
  ######
  # X: a two column matrix containing the coordinate in the plane of points
  ######
  # Output:
  #####
  # Y:points of the border of the smoothed shape described by X

  x_axis <- seq(min(X[, 1]), max(X[, 1]), length.out = prec)
  y_axis <- seq(min(X[, 2]), max(X[, 2]), length.out = prec)

  # Borders
  Y1 <- matrix(X[which(X[, 1] == min(X[, 1]))[1], ], 1, 2)
  for (i in 2:prec) {
    uu <- which(x_axis[i - 1] <= X[, 1] & X[, 1] <= x_axis[i])
    if (length(uu) > 0) {
      new_max <- max(X[uu, 2])
      if (new_max > Y1[nrow(Y1), 2]) {
        Y1 <- rbind(Y1, X[uu[which(X[uu, 2] == new_max)[1]], ])
      }
    }
  }
  Y2 <- matrix(X[which(X[, 1] == max(X[, 1]))[1], ], 1, 2)
  for (i in prec:2) {
    uu <- which(x_axis[i - 1] <= X[, 1] & X[, 1] <= x_axis[i])
    if (length(uu) > 0) {
      new_min <- min(X[uu, 2])
      if (new_min < Y2[nrow(Y2), 2]) {
        Y2 <- rbind(Y2, X[uu[which(X[uu, 2] == new_min)[1]], ])
      }
    }
  }

  # Smoothing
  Y1 <- smooth.spline(Y1, spar = spar)
  Y2 <- smooth.spline(Y2, spar = spar)
  return(cbind(
    c(Y1$x, rev(Y2$x)),
    c(Y1$y, rev(Y2$y))
  ))
}

Bm_MM <- function(G) { # Estimates 2x2 block matrix using method of moment from a graph sample
  #######
  # Input:
  #######
  # - G : a list of adjacency matrices of simple graphs
  ########
  # Output:
  ########
  # - B : The estimated blockmatrix
  #######

  # Estimating densities in oberserved graph
  mu_F <- matrix(NA, length(G), 3)
  for (i in 1:length(G)) {
    A <- G[[i]]
    A2 <- A %*% A
    A3 <- A %*% A2
    mu_F[i, 1] <- sum(A) / (2 * choose(nrow(A), 2)) # Edge density
    mu_F[i, 2] <- sum(diag(A3)) / (6 * choose(nrow(A), 3)) # Triangle density
    mu_F[i, 3] <- (sum(A3 * A) - 4 * sum(choose(diag(A2), 2)) - sum(A)) / (prod(nrow(A):(nrow(A) - 3))) # Square density
  }
  hat_mu_F <- colMeans(mu_F) # Estimator

  # Building objective function function
  obj_funct_1 <- function(theta) {
    return(sum((c(sum(theta^3), sum(theta^4)) - hat_mu_F[-1])^2))
  }
  init <- eigen(apply(array(unlist(G), c(dim(G[[1]]), length(G))), 1:2, mean))$values[1:2] / ncol(G[[1]])
  rez_1 <- optimx(init, obj_funct_1, lower = 0, upper = 1)[1:2]
  hat_a <- sum(rez_1)
  #
  obj_funct_2 <- function(p) {
    return((prod(rez_1) - p * (1 - p) * (hat_a^2 - ((hat_mu_F[1] - (p^2 + (1 - p)^2) * hat_a) / (2 * p * (1 - p)))^2))^2)
  }
  N <- 1e4
  hat_pi <- which.min(obj_funct_2(seq(0, 1, length.out = N))) / N
  hat_b <- (hat_mu_F[1] - (hat_pi^2 + (1 - hat_pi)^2) * hat_a) / (2 * hat_pi * (1 - hat_pi))

  # Output
  #######
  return(list(
    hat_pi = c(hat_pi, 1 - hat_pi),
    hat_B = matrix(c(hat_a, hat_b, hat_b, hat_a), 2, 2)
  ))
}

smpl <- function(A, p1, p2) {
  pmax(
    get.adjacency(erdos.renyi.game(nrow(A), p1), sparse = FALSE) * (1 - A),
    get.adjacency(erdos.renyi.game(nrow(A), p2), sparse = FALSE) * A
  )
}

rG <- function(n, B, P) {
  # Sample a planted partition Bm
  x <- (1:n) / (n + 1)
  get.adjacency(sample_sbm(n, B, c(sum(x < P[1]), sum(x > P[1]))), sparse = FALSE)
}
