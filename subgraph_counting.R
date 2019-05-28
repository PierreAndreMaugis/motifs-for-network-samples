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

require(igraph)
require(Matrix)
ffct <- function(x, k) {
  sapply(x, function(y) {
    prod(y:(y - k + 1))
  })
} # Falling factorial

get.backtrack <- function(A) {
  # Input: The adjacency matrix of a simple graph g
  # Output: the non-backtrakcing matrix of g

  # e is a list of all edges, ordered lexicographicaly (as g is simple)
  e <- get.edgelist(graph.adjacency(A, "upper", diag = FALSE))
  e <- e[do.call(order, lapply(1:ncol(e), function(i) e[, i])), ]
  # Creates a directed edge-list where if ij is an undirected edge in e, then ji suceeds ij in e_dir
  e_dir <- matrix(0, 2 * nrow(e), 2)
  e_dir[seq(1, 2 * nrow(e) - 1, by = 2), ] <- e
  e_dir[seq(2, 2 * nrow(e), by = 2), ] <- e[, c(2, 1)]
  # e_ord is e_dir ordered lexicographicaly. The permutation is saved for future use.
  ord <- do.call(order, lapply(1:ncol(e_dir), function(i) e_dir[, i]))
  e_ord <- e_dir[ord, ]
  # Compute the degree d of each node
  d <- apply(matrix(1:max(e), 1), 2, function(x) {
    length(which(e == x))
  })
  dind <- c(0, cumsum(d))

  # #Creates output matrixIf the network is small, sticking to dense matrices
  if (nrow(e) < 500) {
    b <- matrix(0, nrow(e_dir), nrow(e_dir))
    b_dir <- matrix(0, nrow(e_dir), nrow(e_dir))
  }
  else {
    b <- Matrix(0, nrow(e_dir), nrow(e_dir))
    b_dir <- Matrix(0, nrow(e_dir), nrow(e_dir))
  }
  # Loop over each edge and uses the degree of the tail node to efficiently imput entries in b
  for (i in 1:nrow(e_ord)) {
    j <- e_ord[i, 2] # Label of the tail node
    c <- dind[j] + (1:d[j]) # Nodes connected to j
    cnb <- c[e_ord[c, 2] != e_ord[i, 1]] # Nodes connected to j that are not i
    if (length(cnb) != 0) { # Input matrix b
      b[i, cnb] <- as.double(1.0)
    }
  }
  b_dir[ord, ord] <- b # Reorders the rows and column of b to be as as that of e_dir
  return(b_dir)
}

X_Fs <- function(A) {
  # Input: The edgelist of a simple graph g
  # Output: The number F1 to F15 in g (see list below)
  #
  # Functions to plot the graphs considered using the package igraph
  # F1=graph.ring(3);
  # F2=graph_from_edgelist(t(matrix(c(1,2,1,3,2,3,3,4,4,2),nrow=2)),F);
  # F3=graph_from_edgelist(t(matrix(c(1,2,1,3,2,3,1,4,1,5,4,5),nrow=2)),F);
  # F4=graph_from_edgelist(t(matrix(c(1,2,1,3,2,3,3,4,4,5,5,1),nrow=2)),F);
  # F5=graph_from_edgelist(t(matrix(c(1,2,1,3,3,2,1,4,4,2,1,5,5,2),nrow=2)),F);
  # F6=graph_from_edgelist(t(matrix(c(1,2,1,3,2,3,1,4,4,5,5,6,6,1),nrow=2)),F);
  # F7=graph_from_edgelist(t(matrix(c(1,2,2,3,3,4,4,1,1,3,2,4),nrow=2)),F);
  # F8=graph.ring(4);
  # F9=graph_from_edgelist(t(matrix(c(1,2,1,3,2,3,3,4,4,5,5,1,1,4),nrow=2)),F);
  # F10=graph_from_edgelist(t(matrix(c(1,3,3,2,1,4,4,2,1,5,5,2),nrow=2)),F);
  # F11=graph_from_edgelist(t(matrix(c(1,2,2,3,3,4,4,1,1,5,5,6,6,4),nrow=2)),F);
  # F12=graph_from_edgelist(t(matrix(c(1,2,2,3,3,6,6,1,2,4,4,6,2,5,5,6),nrow=2)),F);
  # F13=graph_from_edgelist(t(matrix(c(1,2,2,3,3,6,6,2,2,5,5,4,4,3,3,1),nrow=2)),F);
  # F14=graph_from_edgelist(t(matrix(c(1,2,2,3,3,4,4,5,5,1,2,6,6,1),nrow=2)),F);
  # F15=graph_from_edgelist(t(matrix(c(1,2,2,3,3,4,4,1,1,5,5,6,6,7,7,1),nrow=2)),F);
  if (sum(A) == 0) {
    return(rep(0, 15))
  }
  b <- get.backtrack(A) * 1.0 # Non-backtracking matrix of g with row and column in specific ordering
  # Location of edges and their reverse
  s1 <- seq(1, nrow(b) - 1, by = 2) # Edges in first direction
  s2 <- seq(2, nrow(b), by = 2) # Edges in reverse direction
  # Index to flip direction of all edges
  ss <- 1:nrow(b)
  ss[s1] <- s2
  ss[s2] <- s1
  # Usefull matrices (we optimized computation by storing a matrix if it is used more than once)
  b2 <- b %*% b
  b3 <- b %*% b2
  b4 <- b %*% b3
  b5 <- b %*% b4
  b_sq <- b2 * t(b2) # b_sq[e,e']=1 if there exists e_1 and e_2 such that ee_1e'e_2 is a closed non-backtracking walk of length 4 (i.e. induce a square)
  b_tr_tr <- diag(b3) %*% t(diag(b3)) # b_tr_tr[e,e'] is the number of pairs of triangles such that the first has e as an edge and the other e' as an edge
  b_sq_tr <- diag(b3) %*% t(diag(b4)) # b_sq_tr[e,e'] is the number of pairs of triangle and squares such that the first has e as an edge and the other e' as an edge
  b_sq_r <- b_sq * b2[, ss] # b_sq_r[e,e']=1 if e and e' are connected by a square with a diagonal starting at the heads of e and e' (see F2 below)
  b_sq_tt <- b_sq * b3[, ss] # b_sq_tt[e,e'] is the number of square with a diagonal of length two (see F10 below) conecting e and e'

  X_F1 <- sum(diag(b3)) / 6
  X_F2 <- sum(Schoose2(diag(b3))) / 2
  X_F3 <- (sum(b * b_tr_tr) - 12 * X_F2 - 6 * X_F1) / 8
  X_F4 <- (sum(diag(b3) * diag(b4)) - 8 * X_F2) / 2
  X_F5 <- sum(Schoose3(diag(b3))) / 2
  X_F6 <- (sum(b * b_sq_tr) - 6 * X_F4 - 24 * X_F5 - 16 * X_F2) / 4
  X_F7 <- sum(b2[s1, s1] * b2[s1, s2] * b2[s2, s1] * b2[s2, s2]) / 6
  X_F8 <- sum(b_sq) / 8
  X_F9 <- sum(b_sq_r * b3) / 2
  X_F10 <- sum(b_sq_tt) / 12
  X_F11 <- (sum(choose(diag(b4), 2)) - 2 * X_F9 - 12 * X_F10 - 12 * X_F7) / 2
  X_F12 <- sum(b_sq * Schoose2(b3[, ss])) / 24
  X_F13 <- sum(b_sq * Schoose2(b3)) / 2
  X_F14 <- (sum(diag(b3) * diag(b5)) - 4 * X_F4 - 2 * X_F9) / 2
  X_F15 <- (sum(b * (diag(b4) %*% t(diag(b4)))) - 8 * X_F8 - 48 * X_F10 - 48 * X_F12 - 12 * X_F11 - 16 * X_F13 - 20 * X_F9 - 72 * X_F7) / 8
  #
  n <- nrow(A)
  return(c(
    X_F1 / (ffct(n, 3) / 6),
    X_F2 / (ffct(n, 4) / 4),
    X_F3 / (ffct(n, 5) / 8),
    X_F4 / (ffct(n, 5) / 2),
    X_F5 / (ffct(n, 5) / 12),
    X_F6 / (ffct(n, 6) / 4),
    X_F7 / (ffct(n, 4) / 24),
    X_F8 / (ffct(n, 4) / 8),
    X_F9 / (ffct(n, 5) / 2),
    X_F10 / (ffct(n, 5) / 12),
    X_F11 / (ffct(n, 6) / 4),
    X_F12 / (ffct(n, 6) / 48),
    X_F13 / (ffct(n, 6) / 4),
    X_F14 / (ffct(n, 6) / 2),
    X_F15 / (ffct(n, 7) / 8)
  ))
}
