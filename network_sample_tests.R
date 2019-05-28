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
####################################
# Graph notations and visulatisation (and basic checks):
####################################
# require(igraph)
# aut<-function(H){count_subgraph_isomorphisms(H,H)};
# X_F<-function(H,g){count_subgraph_isomorphisms(H,g)/aut(H)};
# G_m<-function(...){graph_from_edgelist(t(matrix(c(...),nrow=2)),F)};
#
# #Triangle Variance
# ##################
# tg0=G_m(1,2,2,3,3,1);
# tg1=G_m(1,2,2,3,3,1,1,4,4,5,5,1);
# tg2=G_m(1,2,2,3,3,4,4,1,1,3);
# tg3=G_m(1,2,2,3,3,1,4,5,5,6,6,4);
# par(mfrow=c(1,4));for (k in 0:3){
# eval(parse(text=paste('plot(tg',k,',main=\'tg',k,'\')',sep='')))}
# #
# g = sample_gnp(6,.5);
# u = X_F(tg0,g)^2;
# v = c(X_F(tg0,g),X_F(tg1,g),X_F(tg2,g),X_F(tg3,g));
# sum(v*c(1,2,2,2))-u;
# c(aut(tg0),aut(tg1),aut(tg2));
# #
# #Square Variance
# ################
# sg0=G_m(1,2,2,3,3,4,4,1);
# sg1=G_m(1,2,2,3,3,4,4,1,1,5,5,6,6,7,7,1);
# sg2=G_m(1,2,2,3,3,4,4,1,1,5,5,3,3,6,6,1);
# sg3=G_m(1,2,2,3,3,4,4,1,1,3,3,5,5,6,6,1);
# sg4=G_m(1,2,2,3,3,4,4,1,4,5,5,6,6,1);
# sg5=G_m(1,2,2,3,3,4,4,1,2,4,4,5,5,1);
# sg6=G_m(1,2,2,3,3,4,4,1,1,5,5,3);
# sg7=G_m(1,2,2,3,3,4,4,1,1,3,2,4);
# sg8=G_m(1,2,2,3,3,4,4,1,5,6,6,7,7,8,8,5);
# par(mfrow=c(3,3));for (k in 0:8){
# eval(parse(text=paste('plot(sg',k,',main=\'sg',k,'\')',sep='')))}
# #
# g = sample_gnp(8,.5);
# u = X_F(sg0,g)^2;
# v = array(NA,9);
# v[1] = X_F(sg0,g);v[2] = X_F(sg1,g);v[3] = X_F(sg2,g);
# v[4] = X_F(sg3,g);v[5] = X_F(sg4,g);v[6] = X_F(sg5,g);
# v[7] = X_F(sg6,g);v[8] = X_F(sg7,g);v[9] = X_F(sg8,g);
# sum(v*c(1,2,6,2,2,2,6,6,2))-u;
# c(aut(sg0),aut(sg1),aut(sg2),aut(sg3),aut(sg4),aut(sg5),aut(sg6),aut(sg7));
# #
# #
# #Edge-Triangle Covariance
# #########################
# etg0=G_m(1,2);
# etg1=G_m(1,2,2,3,3,1);
# etg2=G_m(1,2,2,3,3,1,1,4);
# etg3=G_m(1,2,2,3,3,1,4,5);
# par(mfrow=c(1,4));for (k in 0:3){
# eval(parse(text=paste('plot(etg',k,',main=\'etg',k,'\')',sep='')))}
# #
# g = sample_gnp(8,.5);
# u = X_F(etg0,g)*X_F(etg1,g);
# v = c(X_F(etg1,g),X_F(etg2,g),X_F(etg3,g));
# sum(v*c(3,1,1))-u;
# c(aut(etg1),aut(etg2));
# #
# #Edge-Square Covariance
# #######################
# esg0=G_m(1,2);
# esg1=G_m(1,2,2,3,3,4,4,1);
# esg2=G_m(1,2,2,3,3,4,4,1,1,5);
# esg3=G_m(1,2,2,3,3,4,4,1,1,3);
# esg4=G_m(1,2,2,3,3,4,4,1,5,6);
# par(mfrow=c(1,5));for (k in 0:4){
# eval(parse(text=paste('plot(esg',k,',main=\'esg',k,'\')',sep='')))}
# #
# g = sample_gnp(8,.5);
# u = X_F(esg0,g)*X_F(esg1,g);
# v = c(X_F(esg1,g),X_F(esg2,g),X_F(esg3,g),X_F(esg4,g));
# sum(v*c(4,1,1,1))-u;
# c(aut(esg1),aut(esg2),aut(esg3));
# #
# #Triangle-Square Covariance
# ###########################
# tsg0=G_m(1,2,2,3,3,1);
# tsg1=G_m(1,2,2,3,3,4,4,1);
# tsg2=G_m(1,2,2,3,3,4,4,1,1,5,5,6,6,1);
# tsg3=G_m(1,2,2,3,3,4,4,1,1,3,3,5,5,1);
# tsg4=G_m(1,2,2,3,3,4,4,1,1,5,5,2);
# tsg5=G_m(1,2,2,3,3,4,4,1,1,3);
# tsg6=G_m(1,2,2,3,3,4,4,1,5,6,6,7,7,5);
# par(mfrow=c(2,4));for (k in 0:6){
# eval(parse(text=paste('plot(tsg',k,',main=\'tsg',k,'\')',sep='')))}
# #
# g = sample_gnp(7,.5);
# u = X_F(tsg0,g)*X_F(tsg1,g);
# v = c(X_F(tsg2,g),X_F(tsg3,g),X_F(tsg4,g),X_F(tsg5,g),X_F(tsg6,g));
# sum(v*c(1,3,1,2,1))-u;
# c(aut(tsg2),aut(tsg3),aut(tsg4),aut(tsg5));


require(mixtools)
require(SDMTools)
require(igraph)
require(Matrix)
SBm_CR <- function(G, B, P, level = .95, do.plot = TRUE) {
  #######
  # Input:
  #######
  # - G    : a list of adjacency matrices of simple graphs
  # - P    : the probability of being in each group
  # - B    : a |P|x|P| block matrix for the FSBm
  # - level: confidence level for the test
  ########
  # Output:
  ########
  # - hat_mu : The estimated subgraph density
  # - p.value: The p-value for the test
  #######
  # At the top of this file, we provide the graph notation that we use.
  # There, we also present code to first visulise the graphs, and then to check
  # the constant factors used in the formulas for the covariances.

  # Counting C_2,C_3,C_4 in G to estimate hat_mu_F
  ###############################################
  mu_F0 <- matrix(NA, length(G), 3)
  for (i in 1:length(G)) {
    A <- G[[i]]
    A2 <- A %*% A
    A3 <- A %*% A2
    mu_F0[i, 1] <- sum(A) / (2 * choose(nrow(A), 2)) # edge count
    mu_F0[i, 2] <- sum(diag(A3)) / (6 * choose(nrow(A), 3)) # triangle count
    mu_F0[i, 3] <- (sum(A3 * A) - 4 * sum(choose(diag(A2), 2)) - sum(A)) / (prod(nrow(A):(nrow(A) - 3))) # squares count
  }
  hat_mu_F <- colMeans(mu_F0) # Estimator
  #
  # Building mean and covariance under the null
  ############################################
  k <- nrow(B) # Number of blocks
  # Edge, triangle and square densities
  ed <- 0
  tr <- 0
  sq <- 0
  # To store components of the covariance (i.e., densities of graphs that can be built form overlapping copies of
  # edges, triangles and squares)
  edg <- array(0, 2) # Subgraphs that can be built by overlapping an edge and an edge
  tri <- array(0, 3) # Subgraphs that can be built by overlapping a triangle and a triangle
  sqr <- array(0, 8) # Subgraphs that can be built by overlapping a square and a square
  etr <- array(0, 2) # Subgraphs that can be built by overlapping an edge a triangle
  esq <- array(0, 3) # Subgraphs that can be built by overlapping an edge and a square
  tsq <- array(0, 4) # Subgraphs that can be built by overlapping a triangle and a square
  for (i1 in 1:k) {
    # Big imbricated loops to compute densities of each requiered subgraph; for each vertex in each subgraph, for
    # each possible block allocation, compure probability of said allocation and joint edge occurence. As edge and
    # block allocations are independence, this is a simple product.
    for (i2 in 1:k) {
      ed <- ed + P[i1] * P[i2] * B[i1, i2]
      edg[2] <- edg[2] + P[i1] * P[i2] * B[i1, i2]
      for (i3 in 1:k) {
        tr <- tr + P[i1] * P[i2] * P[i3] * B[i1, i2] * B[i2, i3] * B[i3, i1]
        edg[1] <- edg[1] + P[i1] * P[i2] * P[i3] * B[i1, i2] * B[i2, i3]
        tri[3] <- tri[3] + P[i1] * P[i2] * P[i3] * B[i1, i2] * B[i2, i3] * B[i3, i1]
        etr[2] <- etr[2] + P[i1] * P[i2] * P[i3] * B[i1, i2] * B[i2, i3] * B[i3, i1]
        for (i4 in 1:k) {
          sq <- sq + P[i1] * P[i2] * P[i3] * P[i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1]
          tri[2] <- tri[2] + P[i1] * P[i2] * P[i3] * P[i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i3]
          sqr[7] <- sqr[7] + P[i1] * P[i2] * P[i3] * P[i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i3] * B[i2, i4]
          sqr[8] <- sqr[8] + P[i1] * P[i2] * P[i3] * P[i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1]
          etr[1] <- etr[1] + P[i1] * P[i2] * P[i3] * P[i4] * B[i1, i2] * B[i2, i3] * B[i3, i1] * B[i1, i4]
          esq[2] <- esq[2] + P[i1] * P[i2] * P[i3] * P[i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i3]
          esq[3] <- esq[3] + P[i1] * P[i2] * P[i3] * P[i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1]
          tsq[4] <- tsq[4] + P[i1] * P[i2] * P[i3] * P[i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i3]
          for (i5 in 1:k) {
            tri[1] <- tri[1] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * B[i1, i2] * B[i2, i3] * B[i3, i1] * B[i1, i4] * B[i4, i5] *
              B[i5, i1]
            sqr[5] <- sqr[5] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i5] *
              B[i5, i4] * B[i4, i2]
            sqr[6] <- sqr[6] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i5] *
              B[i5, i3]
            esq[1] <- esq[1] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i5]
            tsq[2] <- tsq[2] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i3] *
              B[i3, i5] * B[i5, i1]
            tsq[3] <- tsq[3] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i5] *
              B[i5, i2]
            for (i6 in 1:k) {
              sqr[2] <- sqr[2] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * P[i6] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] *
                B[i1, i5] * B[i5, i3] * B[i3, i6] * B[i6, i1]
              sqr[3] <- sqr[3] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * P[i6] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] *
                B[i1, i3] * B[i3, i5] * B[i5, i6] * B[i6, i1]
              sqr[4] <- sqr[4] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * P[i6] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] *
                B[i1, i5] * B[i5, i6] * B[i6, i4]
              tsq[1] <- tsq[1] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * P[i6] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] *
                B[i1, i5] * B[i5, i6] * B[i6, i1]
              for (i7 in 1:k) {
                sqr[1] <- sqr[1] + P[i1] * P[i2] * P[i3] * P[i4] * P[i5] * P[i6] * P[i7] * B[i1, i2] * B[i2, i3] * B[i3, i4] *
                  B[i4, i1] * B[i1, i5] * B[i5, i6] * B[i6, i7] * B[i7, i1]
              }
            }
          }
        }
      }
    }
  }
  mu_F <- c(ed, tr, sq)
  # Centring; i.e., removing the term from non-overlapping terms
  edg <- edg - ed^2
  tri <- tri - tr^2
  sqr <- sqr - sq^2
  etr <- etr - ed * tr
  esq <- esq - ed * sq
  tsq <- tsq - tr * sq

  # Utilities
  n <- sapply(G, nrow) # Network sizes
  ffct <- function(x, k) {
    sapply(x, function(y) {
      prod(y:(y - k + 1))
    })
  } # Falling factorial

  # The covariances (see top file comments for graph naming and plots; the factors c_H are computed by hand)
  sgm <- matrix(NA, 3, 3) # Cov matrixfor worse case (ER case; see SI)
  sgm[1, 1] <- sum((2 * (ffct(n, 3) / 2) * edg[1] + # Chery
    1 * choose(n, 2) * edg[2]) / (choose(n, 2)^2)) # Edge
  sgm[2, 2] <- sum((2 * (ffct(n, 5) / 8) * tri[1] + # Butterfly; tg1
    2 * (ffct(n, 4) / 4) * tri[2] + # Diamond; tg2
    choose(n, 3) * tri[3]) / (choose(n, 3)^2)) # Triangle; tg0
  sgm[3, 3] <- sum((2 * (ffct(n, 7) / 8) * sqr[1] + # sg1
    6 * (ffct(n, 6) / 48) * sqr[2] + # sg2
    2 * (ffct(n, 6) / 4) * sqr[3] + # sg3
    2 * (ffct(n, 6) / 4) * sqr[4] + # sg4
    2 * (ffct(n, 5) / 2) * sqr[5] + # sg5
    6 * (ffct(n, 5) / 12) * sqr[6] + # sg6
    6 * (ffct(n, 4) / 24) * sqr[7] + # sg7
    (ffct(n, 4) / 8) * sqr[8]) / (ffct(n, 4) / 8)^2) # Square
  sgm[1, 2] <- sgm[2, 1] <- sum((
    1 * (ffct(n, 4) / 2) * etr[1] + # Shovel
      3 * choose(n, 3) * etr[2]) / (choose(n, 3) * choose(n, 2))) # Triangle
  sgm[1, 3] <- sgm[3, 1] <- sum((
    1 * (ffct(n, 5) / 2) * esq[1] + # Banner
      1 * (ffct(n, 4) / 4) * esq[2] + # Diamond
      4 * (ffct(n, 4) / 8) * esq[3]) / ((ffct(n, 4) / 8) * choose(n, 2))) # Square
  sgm[2, 3] <- sgm[3, 2] <- sum((
    1 * (ffct(n, 6) / 4) * tsq[1] + # tsg2
      3 * (ffct(n, 5) / 12) * tsq[2] + # tsg3
      1 * (ffct(n, 5) / 2) * tsq[3] + # tsg4
      2 * (ffct(n, 4) / 4) * tsq[4]) / ((ffct(n, 4) / 8) * choose(n, 3))) # tsg5

  # Building the ellipse and p-value
  #################################
  conf_region <- list()
  sgm_t <- list()
  for (i in 3:1) {
    j <- (3:1)[i] # Numbering that match the rest of the notation, but not the loop
    sgm_t[[j]] <- sgm[-i, -i] / length(G)^2 # generator of auxillary spheres
    conf_region[[j]] <- ellipse(mu_F[-i], sgm_t[[j]], alpha = 1 - level, draw = FALSE)
  }
  sgm_inv <- solve(sgm)
  stat <- (hat_mu_F - mu_F) %*% sgm_inv %*% (hat_mu_F - mu_F) * length(G)^2
  p.value <- 1 - pchisq(stat, 3)

  # Plotting
  #########
  if (do.plot) {
    par(mfrow = c(1, 3))
    labs <- c("C2", "C3", "C4")
    for (i in 3:1) {
      j <- (3:1)[i]
      plot(NA,
        ylim = pmax(range(conf_region[[j]][, 2], hat_mu_F[-i][2]) * c(.9, 1.1), 0),
        xlim = pmax(range(conf_region[[j]][, 1], hat_mu_F[-i][1]) * c(.9, 1.1), 0),
        xlab = (labs[-i])[1], ylab = (labs[-i])[2]
      )
      points(mu_F0[, -i][, 1],
        mu_F0[, -i][, 2],
        cex = .2,
        pch = 2
      )
      polygon(conf_region[[j]],
        col = rgb(1, 0, 0, .1),
        border = rgb(1, 0, 0, 1),
        lwd = .5
      )
      points(mu_F[-i][1],
        mu_F[-i][2],
        pch = 16,
        col = rgb(1, 0, 0, 1)
      )
      points(hat_mu_F[-i][1],
        hat_mu_F[-i][2],
        pch = 3
      )
      points(
        hat_mu_F[-i][1],
        hat_mu_F[-i][2]
      )
    }
  }

  # Output
  #######
  return(list(
    hat_mu = hat_mu_F,
    p.value = p.value,
    dst = sqrt(sum((hat_mu_F - mu_F)^2)),
    conf_region = conf_region,
    mu_F0 = mu_F0,
    mu_F = mu_F,
    stat = stat,
    sgm = sgm
  ))
}

FSBm_CR0 <- function(G, B, level = .95, by = .01, prec = .01, do.plot = TRUE) {
  #######
  # Input:
  #######
  # - G    : a list of adjacency matrices of simple graphs
  # - B    : a block matrix for the FSBm
  # - level: confidence level for the test
  # - by   : precision of embedding shape and convex hull
  # - prec : precision for p.value
  ########
  # Output:
  ########
  # - hat_mu : The estimated subgraph density
  # - p.value: The p-value for the test
  #######
  # At the top of this file, we provide the graph notation that we use.
  # There, we also present code to first visulise the graphs, and then to check
  # the constant factors used in the formulas for the covariances.

  # Counting C_2,C_3,C_4 in G to estimate hat_mu_F
  ###############################################
  mu_F0 <- matrix(NA, length(G), 3)
  for (i in 1:length(G)) {
    A <- G[[i]]
    A2 <- A %*% A
    A3 <- A %*% A2
    mu_F0[i, 1] <- sum(A) / (2 * choose(nrow(A), 2)) # Edge count
    mu_F0[i, 2] <- sum(diag(A3)) / (6 * choose(nrow(A), 3)) # Triangle count
    mu_F0[i, 3] <- (sum(A3 * A) - 4 * sum(choose(diag(A2), 2)) - sum(A)) / (prod(nrow(A):(nrow(A) - 3))) # Square Count
  }
  hat_mu_F <- colMeans(mu_F0) # Estimator

  # Building paremtrization of FSBM subgraph densities set Pi
  ###################################################
  k <- nrow(B) # Number of blocks
  x <- expand.grid(replicate(k - 1, seq(0, 1, by = by), simplify = FALSE))
  x <- cbind(x, 1 - rowSums(x))
  x <- x[which(x[, k] >= 0), ]
  x <- t(apply(x, 1, sort))
  x <- x[duplicated(x) == 0, ]

  # Building mean and covariance under the null
  ############################################
  # Edge, triangle and square densities
  ed <- array(0, length(x))
  tr <- array(0, length(x))
  sq <- array(0, length(x))
  # To store components of the covariance (i.e., densities of graphs that can be built form overlapping copies of
  # edges, triangles and squares)
  edg <- array(0, c(length(x), 2))
  tri <- array(0, c(length(x), 3))
  sqr <- array(0, c(length(x), 8))
  etr <- array(0, c(length(x), 2))
  esq <- array(0, c(length(x), 3))
  tsq <- array(0, c(length(x), 4))
  #
  for (i1 in 1:k) {
    # Big imbricated loops to compute densities of each requiered subgraph; for each vertex in each subgraph, for
    # each possible block allocation, compure probability of said allocation and joint edge occurence. As edge and
    # block allocations are independence, this is a simple product.
    for (i2 in 1:k) {
      ed <- ed + x[, i1] * x[, i2] * B[i1, i2]
      edg[, 2] <- edg[, 2] + x[, i1] * x[, i2] * B[i1, i2]
      for (i3 in 1:k) {
        tr <- tr + x[, i1] * x[, i2] * x[, i3] * B[i1, i2] * B[i2, i3] * B[i3, i1]
        edg[, 1] <- edg[, 1] + x[, i1] * x[, i2] * x[, i3] * B[i1, i2] * B[i2, i3]
        tri[, 3] <- tri[, 3] + x[, i1] * x[, i2] * x[, i3] * B[i1, i2] * B[i2, i3] * B[i3, i1]
        etr[, 2] <- etr[, 2] + x[, i1] * x[, i2] * x[, i3] * B[i1, i2] * B[i2, i3] * B[i3, i1]
        for (i4 in 1:k) {
          sq <- sq + x[, i1] * x[, i2] * x[, i3] * x[, i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1]
          tri[, 2] <- tri[, 2] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i3]
          sqr[, 7] <- sqr[, 7] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i3] *
            B[i2, i4]
          sqr[, 8] <- sqr[, 8] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1]
          etr[, 1] <- etr[, 1] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * B[i1, i2] * B[i2, i3] * B[i3, i1] * B[i1, i4]
          esq[, 2] <- esq[, 2] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i3]
          esq[, 3] <- esq[, 3] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1]
          tsq[, 4] <- tsq[, 4] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] * B[i1, i3]
          for (i5 in 1:k) {
            tri[, 1] <- tri[, 1] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * B[i1, i2] * B[i2, i3] * B[i3, i1] * B[i1, i4] *
              B[i4, i5] * B[i5, i1]
            sqr[, 5] <- sqr[, 5] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] *
              B[i1, i5] * B[i5, i4] * B[i4, i2]
            sqr[, 6] <- sqr[, 6] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] *
              B[i1, i5] * B[i5, i3]
            esq[, 1] <- esq[, 1] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] *
              B[i1, i5]
            tsq[, 2] <- tsq[, 2] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] *
              B[i1, i3] * B[i3, i5] * B[i5, i1]
            tsq[, 3] <- tsq[, 3] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1] *
              B[i1, i5] * B[i5, i2]
            for (i6 in 1:k) {
              sqr[, 2] <- sqr[, 2] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * x[, i6] * B[i1, i2] * B[i2, i3] * B[i3, i4] *
                B[i4, i1] * B[i1, i5] * B[i5, i3] * B[i3, i6] * B[i6, i1]
              sqr[, 3] <- sqr[, 3] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * x[, i6] * B[i1, i2] * B[i2, i3] * B[i3, i4] *
                B[i4, i1] * B[i1, i3] * B[i3, i5] * B[i5, i6] * B[i6, i1]
              sqr[, 4] <- sqr[, 4] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * x[, i6] * B[i1, i2] * B[i2, i3] * B[i3, i4] *
                B[i4, i1] * B[i1, i5] * B[i5, i6] * B[i6, i4]
              tsq[, 1] <- tsq[, 1] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * x[, i6] * B[i1, i2] * B[i2, i3] * B[i3, i4] *
                B[i4, i1] * B[i1, i5] * B[i5, i6] * B[i6, i1]
              for (i7 in 1:k) {
                sqr[, 1] <- sqr[, 1] + x[, i1] * x[, i2] * x[, i3] * x[, i4] * x[, i5] * x[, i6] * x[, i7] * B[i1, i2] * B[i2, i3] *
                  B[i3, i4] * B[i4, i1] * B[i1, i5] * B[i5, i6] * B[i6, i7] * B[i7, i1]
              }
            }
          }
        }
      }
    }
  }

  # Removing replicates and centering
  super_mu <- cbind(ed, tr, sq, edg, tri, sqr, etr, esq, tsq)
  u <- which(duplicated(super_mu) == 0)
  mu_F <- cbind(ed, tr, sq)[u, ]
  edg <- sweep(edg, 1, ed^2)[u, ]
  tri <- sweep(tri, 1, tr^2)[u, ]
  sqr <- sweep(sqr, 1, sq^2)[u, ]
  etr <- sweep(etr, 1, ed * tr)[u, ]
  esq <- sweep(esq, 1, ed * sq)[u, ]
  tsq <- sweep(tsq, 1, tr * sq)[u, ]

  # Utilities
  n <- sapply(G, nrow) # network sizes
  ffct <- function(x, k) {
    sapply(x, function(y) {
      prod(y:(y - k + 1))
    })
  } # Falling factorial

  # The covariances (see top file comments for graph naming and plots; the factors c_H are computed by hand)
  sgm <- array(NA, c(3, 3, length(u))) # Cov matrices
  for (l in 1:length(u)) {
    sgm[1, 1, l] <- sum((2 * (ffct(n, 3) / 2) * edg[l, 1] + # Chery
      1 * choose(n, 2) * edg[l, 2]) / (choose(n, 2)^2)) # Edge
    sgm[2, 2, l] <- sum((2 * (ffct(n, 5) / 8) * tri[l, 1] + # Butterfly; tg1
      2 * (ffct(n, 4) / 4) * tri[l, 2] + # Diamond; tg2
      choose(n, 3) * tri[l, 3]) / (choose(n, 3)^2)) # Triangle; tg0
    sgm[3, 3, l] <- sum((2 * (ffct(n, 7) / 8) * sqr[l, 1] + # sg1
      6 * (ffct(n, 6) / 48) * sqr[l, 2] + # sg2
      2 * (ffct(n, 6) / 4) * sqr[l, 3] + # sg3
      2 * (ffct(n, 6) / 4) * sqr[l, 4] + # sg4
      2 * (ffct(n, 5) / 2) * sqr[l, 5] + # sg5
      6 * (ffct(n, 5) / 12) * sqr[l, 6] + # sg6
      6 * (ffct(n, 4) / 24) * sqr[l, 7] + # sg7
      (ffct(n, 4) / 8) * sqr[l, 8]) / (ffct(n, 4) / 8)^2) # Square
    sgm[1, 2, l] <- sgm[2, 1, l] <- sum((
      1 * (ffct(n, 4) / 2) * etr[l, 1] + # Shovel
        3 * choose(n, 3) * etr[l, 2]) / (choose(n, 3) * choose(n, 2))) #  Triangle
    sgm[1, 3, l] <- sgm[3, 1, l] <- sum((
      1 * (ffct(n, 5) / 2) * esq[l, 1] + # Banner
        1 * (ffct(n, 4) / 4) * esq[l, 2] + # Diamond
        4 * (ffct(n, 4) / 8) * esq[l, 3]) / ((ffct(n, 4) / 8) * choose(n, 2))) # Square
    sgm[2, 3, l] <- sgm[3, 2, l] <- sum((
      1 * (ffct(n, 6) / 4) * tsq[l, 1] + # tsg2
        3 * (ffct(n, 5) / 12) * tsq[l, 2] + # tsg3
        1 * (ffct(n, 5) / 2) * tsq[l, 3] + # tsg4
        2 * (ffct(n, 4) / 4) * tsq[l, 4]) / ((ffct(n, 4) / 8) * choose(n, 3))) # tsg5
  }

  # Building the ellipses and p-values
  ###################################
  conf_region <- list()
  sgm_t <- list()
  for (i in 3:1) {
    j <- (3:1)[i] # Numbering that match the rest of the notation, but not the loop
    conf_region[[j]] <- c(NA, NA)
    for (l in 1:length(u)) {
      sgm_t[[j]] <- sgm[-i, -i, l] / length(G)^2 # Generator of ellipse
      conf_region[[j]] <- rbind(
        conf_region[[j]],
        ellipse(mu_F[l, -i],
          sgm_t[[j]],
          alpha = 1 - level,
          npoints = 50,
          draw = FALSE
        )
      )
    }
  }
  p.values <- array(NA, length(u))
  for (l in 1:length(u)) {
    sgm_inv <- solve(sgm[, , l])
    p.values[l] <- 1 - pchisq((hat_mu_F - mu_F[l, ]) %*% sgm_inv %*% (hat_mu_F - mu_F[l, ]) * length(G)^2, 3)
  }

  # Plotting
  #########
  if (do.plot) {
    par(mfrow = c(1, 3))
    labs <- c("C2", "C3", "C4")
    for (i in 3:1) {
      j <- (3:1)[i]
      plot(NA,
        ylim = pmax(
          range(conf_region[[j]][, 2],
            mu_F[, -i][, 2],
            hat_mu_F[-i][2],
            na.rm = TRUE
          ) * c(.9, 1.1),
          0
        ),
        xlim = pmax(
          range(conf_region[[j]][, 1],
            mu_F[, -i][, 1],
            hat_mu_F[-i][1],
            na.rm = TRUE
          ) * c(.9, 1.1),
          0
        ),
        xlab = (labs[-i])[1],
        ylab = (labs[-i])[2],
        axes = FALSE
      )
      for (l in 1:length(u)) {
        sgm_t[[j]] <- sgm[-i, -i, l] / length(G)^2
        ellps <- ellipse(mu_F[l, -i],
          sgm_t[[j]],
          alpha = 1 - level,
          draw = FALSE
        )
        polygon(ellps, col = rgb(1, .5, .5), border = NA)
      }
      polygon(emb_contour(mu_F[, -i]), col = 2, border = 2, lwd = 2)
      polygon(matrix(c(0, 0, -1, -1, -1, 1, 1, -1), 4, 2), col = "white", border = NA)
      polygon(matrix(c(-1, -1, 1, 1, 0, -1, -1, 0), 4, 2), col = "white", border = NA)
      points(hat_mu_F[-i][1], hat_mu_F[-i][2], pch = 3)
      points(hat_mu_F[-i][1], hat_mu_F[-i][2])
      axis(1)
      axis(2)
      box()
    }
  }

  # Output
  #######
  return(list(
    hat_mu = hat_mu_F,
    emb_shape = mu_F,
    p.value = max(p.values),
    conf_region = conf_region,
    mu_F0 = mu_F0,
    mu_F = mu_F,
    sgm = sgm,
    u = u
  ))
}

FSBm_CR <- function(G, B, level = .95, by = 1 / (2 * nrow(B)), prec = .01, do.plot = TRUE) {
  #######
  # Input:
  #######
  # - G    : a list of adjacency matrices of simple graphs
  # - B    : a block matrix for the FSBm
  # - level: confidence level for the test
  # - by   : precision of embedding shape and convex hull
  # - prec : precision for p.value
  ########
  # Output:
  ########
  # - hat_mu : The estimated subgraph density
  # - p.value: The p-value for the test
  #######
  # At the end of the function, we provide the graph notation that we use.
  # There, we also present code to first visulise the graphs, and then to check
  # the constant factors used in the formulas for the covariances.
  #######

  # Counting C_2,C_3,C_4 in G to estimate hat_mu_F
  ###############################################
  mu_F <- matrix(NA, length(G), 3)
  for (i in 1:length(G)) {
    A <- G[[i]]
    A2 <- A %*% A
    A3 <- A %*% A2
    mu_F[i, 1] <- sum(A) / (2 * choose(nrow(A), 2)) # Edge count
    mu_F[i, 2] <- sum(diag(A3)) / (6 * choose(nrow(A), 3)) # Triangle count
    mu_F[i, 3] <- (sum(A3 * A) - 4 * sum(choose(diag(A2), 2)) - sum(A)) / (prod(nrow(A):(nrow(A) - 3))) # Square Count
  }
  hat_mu_F <- colMeans(mu_F) # Estimator

  # Building embedding shapes
  ##########################
  k <- nrow(B) # Number of blocks
  x <- expand.grid(replicate(k - 1, seq(0, 1, by = by), simplify = FALSE)) # Paremtrization of Pi
  x <- cbind(x, 1 - rowSums(x))
  x <- x[which(x[, k] >= 0), ]
  x <- t(apply(x, 1, sort))
  x <- x[duplicated(x) == 0, ]
  ed <- 0
  tr <- 0
  sq <- 0
  for (i1 in 1:k) {
    for (i2 in 1:k) {
      ed <- ed + x[, i1] * x[, i2] * B[i1, i2]
      for (i3 in 1:k) {
        tr <- tr + x[, i1] * x[, i2] * x[, i3] * B[i1, i2] * B[i2, i3] * B[i3, i1]
        for (i4 in 1:k) {
          sq <- sq + x[, i1] * x[, i2] * x[, i3] * x[, i4] * B[i1, i2] * B[i2, i3] * B[i3, i4] * B[i4, i1]
        }
      }
    }
  }
  emd_shape_1 <- cbind(ed, tr) # Embedding shape edge triangles
  emd_shape_2 <- cbind(ed, sq) # Embedding shape edge squares
  emd_shape_3 <- cbind(tr, sq) # Embedding shape triangles squares
  emd_shape <- list(emd_shape_1, emd_shape_2, emd_shape_3)

  # Building convex hulls
  ######################
  ch1 <- emd_shape_1[chull(emd_shape_1), ] # Convex hull
  ch2 <- emd_shape_2[chull(emd_shape_2), ] # Convex hull
  ch3 <- emd_shape_3[chull(emd_shape_3), ] # Convex hull
  ch0 <- list(ch1, ch2, ch3) # For convenient access

  # Building auxillary spheres for each vertex of the hulls
  ########################################################
  p <- max(diag(B)) # Denser block
  n <- sapply(G, nrow) # Network sizes
  sgm <- matrix(NA, 3, 3) # Cov matrix
  ffct <- function(x, k) {
    sapply(x, function(y) {
      prod(y:(y - k + 1))
    })
  } # Falling factorial
  # The covariances (see top file comments for graph naming and plots; the factors c_H are computed by hand)
  sgm[1, 1] <- sum((2 * (ffct(n, 3) / 2) * (p^2) + # Chery
    1 * choose(n, 2) * (p)) / (choose(n, 2)^2)) # Edge
  sgm[2, 2] <- sum((2 * (ffct(n, 5) / 8) * (p^6) + # Butterfly; tg1
    2 * (ffct(n, 4) / 4) * (p^5) + # Diamond; tg2
    choose(n, 3) * (p^3)) / (choose(n, 3)^2)) # Triangle; tg0
  sgm[3, 3] <- sum((2 * (ffct(n, 7) / 8) * (p^8) + # sg1 (see end comments)
    6 * (ffct(n, 6) / 48) * (p^8) + # sg2
    2 * (ffct(n, 6) / 4) * (p^8) + # sg3
    2 * (ffct(n, 6) / 4) * (p^7) + # sg4
    2 * (ffct(n, 5) / 2) * (p^7) + # sg5
    6 * (ffct(n, 5) / 12) * (p^6) + # sg6
    6 * (ffct(n, 4) / 24) * (p^6) + # sg7
    (ffct(n, 4) / 8) * (p^4)) / (ffct(n, 4) / 8)^2) # Square
  sgm[1, 2] <- sgm[2, 1] <- sum((
    1 * (ffct(n, 4) / 2) * (p^4) + # Shovel
      3 * choose(n, 3) * (p^3)) / (choose(n, 3) * choose(n, 2))) # Triangle
  sgm[1, 3] <- sgm[3, 1] <- sum((
    1 * (ffct(n, 5) / 2) * (p^5) + # Banner
      1 * (ffct(n, 4) / 4) * (p^5) + # Diamond
      4 * (ffct(n, 4) / 8) * (p^4)) / ((ffct(n, 4) / 8) * choose(n, 2))) # Square
  sgm[2, 3] <- sgm[3, 2] <- sum((
    1 * (ffct(n, 6) / 4) * (p^7) + # tsg2
      3 * (ffct(n, 5) / 12) * (p^7) + # tsg3
      1 * (ffct(n, 5) / 2) * (p^6) + # tsg4
      2 * (ffct(n, 4) / 4) * (p^5)) / ((ffct(n, 4) / 8) * choose(n, 3))) # tsg5

  # Building the ellipses
  aux_sph <- list()
  sgm_t <- list()
  for (i in 3:1) {
    j <- (3:1)[i] # Numbering that match the rest of the notation, but not the loop
    l <- max(eigen(sgm[-i, -i])$values) / length(G)^2 # Radius of the auxillary spheres
    sgm_t[[j]] <- matrix(c(l, 0, 0, l), 2, 2) # Generator of auxillary spheres
    aux_sph0 <- c() # Temporary holder for auxillary spheres
    for (y in 1:nrow(ch0[[j]])) {
      aux_sph0 <- rbind(
        aux_sph0,
        ellipse((ch0[[j]])[y, ],
          sgm_t[[j]],
          alpha = 1 - level,
          draw = FALSE
        )
      )
    }
    aux_sph[[j]] <- aux_sph0[chull(aux_sph0), ] # Convex hull of auxillary spheres
  }
  aux_sph1 <- aux_sph[[1]] # For edge-triangle
  aux_sph2 <- aux_sph[[2]] # For edge-square
  aux_sph3 <- aux_sph[[3]] # For triangle-square

  # Computing the p-value (inefficient loop)
  #########################################
  p.val <- 1 - prec
  stp <- FALSE
  while (stp == FALSE && p.val > 0) {
    sph_pval <- list(c(), c(), c()) # Auxillary spheres for p.value
    for (i in 3:1) {
      j <- (3:1)[i]
      for (y in 1:nrow(ch0[[j]])) {
        sph_pval[[j]] <- rbind(
          sph_pval[[j]],
          ellipse((ch0[[j]])[y, ],
            sgm_t[[j]],
            alpha = 1 - p.val,
            draw = FALSE
          )
        )
      }
      sph_pval[[j]] <- sph_pval[[j]][chull(sph_pval[[j]]), ] # Convex hull of current spheres
      if (pnt.in.poly(matrix(hat_mu_F[-i], 1, 2), sph_pval[[j]])$pip == 0) {
        stp <- TRUE
      }
    }
    if (stp == FALSE) {
      p.val <- p.val - prec
    }
  }

  # Plotting
  #########
  if (do.plot) {
    par(mfrow = c(1, 3))
    labs <- c("C2", "C3", "C4")
    for (i in 3:1) {
      j <- (3:1)[i]
      plot(NA,
        ylim = pmax(range(
          aux_sph[[j]][, 2],
          hat_mu_F[-i][2]
        ) * c(.9, 1.1), 0),
        xlim = pmax(range(
          aux_sph[[j]][, 1],
          hat_mu_F[-i][1]
        ) * c(.9, 1.1), 0),
        xlab = (labs[-i])[1],
        ylab = (labs[-i])[2]
      )
      # points(mu_F[,-i][,1],mu_F[,-i][,2],cex=.2,pch=2)
      polygon(aux_sph[[j]],
        col = rgb(1, 0, 0, .1),
        border = rgb(1, 0, 0, 1),
        lwd = .5
      )
      polygon(ch0[[j]],
        col = rgb(1, 0, 0, .3),
        border = 1
      )
      # lines(emd_shape,type='l',lwd=2,col=rgb(1,0,0,1))
      polygon(ch0[[j]],
        border = 1,
        lwd = 2
      )
      points(hat_mu_F[-i][1],
        hat_mu_F[-i][2],
        pch = 3
      )
      points(
        hat_mu_F[-i][1],
        hat_mu_F[-i][2]
      )
    }
  }

  # Output
  #######
  return(list(
    hat_mu = hat_mu_F,
    p.value = 1 - max(0, p.val),
    mu_F = mu_F,
    aux_sph = aux_sph,
    emd_shape = emd_shape,
    ch0 = ch0
  ))
}

TSMPL_CR <- function(G1, G2) {
  #######
  # Input:
  #######
  # - G1,G2: a list of adjacency matrices of simple graphs
  ########
  # Output:
  ########
  # - hat_mu : The difference in estimated subgraph densities
  # - hat_sgm: The estimated covariance matrix
  # - p.value: The p-value for the test

  #######
  # Counting rooted C_2,C_3,C_4 in G to estimate means
  ###################################################
  hat_mu_H <- matrix(NA, length(G1) + length(G2), 15 + 4)
  mu_F1 <- matrix(NA, length(G1), 3)
  mu_F2 <- matrix(NA, length(G2), 3)
  for (i in 1:length(G1)) {
    A <- G1[[i]]
    A2 <- A %*% A
    A3 <- A %*% A2
    A4 <- A %*% A3
    n <- nrow(A)
    edgs <- diag(A2)
    tris <- diag(A3) / 2
    sqrs <- (diag(A4) - rowSums(A2) - diag(A2) * pmax(0, (diag(A2) - 1))) / 2
    #
    mu_F1[i, 1] <- sum(edgs / 2) / (ffct(n, 2) / 2) # Edge count
    mu_F1[i, 2] <- sum(tris / 3) / (ffct(n, 3) / 6) # Triangle count
    mu_F1[i, 3] <- sum(sqrs / 4) / (ffct(n, 4) / 8) # Squares count
    #
    hat_mu_H[i, 1:15] <- X_Fs(A)
    hat_mu_H[i, 16] <- sum(choose(edgs, 2)) / (ffct(n, 3) / 2)
    hat_mu_H[i, 17] <- sum(edgs / 2) / (ffct(n, 2) / 2)
    hat_mu_H[i, 18] <- sum(pmax(edgs - 2, 0) * tris) / (ffct(n, 4) / 2)
    hat_mu_H[i, 19] <- sum(pmax(edgs - 2, 0) * sqrs) / (ffct(n, 5) / 2)
  }
  for (i in 1:length(G2)) {
    A <- G2[[i]]
    A2 <- A %*% A
    A3 <- A %*% A2
    A4 <- A %*% A3
    n <- nrow(A)
    edgs <- diag(A2)
    tris <- diag(A3) / 2
    sqrs <- (diag(A4) - rowSums(A2) - diag(A2) * pmax(0, (diag(A2) - 1))) / 2
    #
    mu_F2[i, 1] <- sum(edgs / 2) / (ffct(n, 2) / 2) # Edge count
    mu_F2[i, 2] <- sum(tris / 3) / (ffct(n, 3) / 6) # Triangle count
    mu_F2[i, 3] <- sum(sqrs / 4) / (ffct(n, 4) / 8) # Squares count
    #
    hat_mu_H[length(G1) + i, 1:15] <- X_Fs(A)
    hat_mu_H[length(G1) + i, 16] <- sum(choose(edgs, 2)) / (ffct(n, 3) / 2)
    hat_mu_H[length(G1) + i, 17] <- sum(edgs / 2) / (ffct(n, 2) / 2)
    hat_mu_H[length(G1) + i, 18] <- sum(pmax(edgs - 2, 0) * tris) / (ffct(n, 4) / 2)
    hat_mu_H[length(G1) + i, 19] <- sum(pmax(edgs - 2, 0) * sqrs) / (ffct(n, 5) / 2)
  }
  mu_H <- colMeans(hat_mu_H)
  mu_F <- colMeans(rbind(mu_F1, mu_F2))
  #
  n1 <- sapply(G1, nrow)
  n2 <- sapply(G2, nrow) # Network sizes
  sgm1 <- sgm2 <- matrix(NA, 3, 3) # Cov matrix
  # The covariances (see top file comments for graph naming and plots; the factors c_H are computed by hand)
  sgm1[1, 1] <- sum((2 * (ffct(n1, 3) / 2) * (mu_H[16] - mu_F[1]^2) + # Chery
    1 * choose(n1, 2) * (mu_H[17] - mu_F[1]^2)) / (choose(n1, 2)^2)) # Edge
  sgm1[2, 2] <- sum((2 * (ffct(n1, 5) / 8) * (mu_H[3] - mu_F[2]^2) + # Butterfly; tg1
    2 * (ffct(n1, 4) / 4) * (mu_H[2] - mu_F[2]^2) + # Diamond; tg2
    choose(n1, 3) * (mu_H[1] - mu_F[2]^2)) / (choose(n1, 3)^2)) # Triangle; tg0
  sgm1[3, 3] <- sum((2 * (ffct(n1, 7) / 8) * (mu_H[15] - mu_F[3]^2) + # sg1 (see end comments)
    6 * (ffct(n1, 6) / 48) * (mu_H[12] - mu_F[3]^2) + # sg2
    2 * (ffct(n1, 6) / 4) * (mu_H[13] - mu_F[3]^2) + # sg3
    2 * (ffct(n1, 6) / 4) * (mu_H[11] - mu_F[3]^2) + # sg4
    2 * (ffct(n1, 5) / 2) * (mu_H[9] - mu_F[3]^2) + # sg5
    6 * (ffct(n1, 5) / 12) * (mu_H[10] - mu_F[3]^2) + # sg6
    6 * (ffct(n1, 4) / 24) * (mu_H[7] - mu_F[3]^2) + # sg7
    1 * (ffct(n1, 4) / 8) * (mu_H[8] - mu_F[3]^2)) / (ffct(n1, 4) / 8)^2) # Square
  sgm1[1, 2] <- sgm1[2, 1] <- sum((
    1 * (ffct(n1, 4) / 2) * (mu_H[18] - mu_F[1] * mu_F[2]) + # Shovel
      3 * choose(n1, 3) * (mu_H[1] - mu_F[1] * mu_F[2])) / (choose(n1, 3) * choose(n1, 2))) # Triangle
  sgm1[1, 3] <- sgm1[3, 1] <- sum((
    1 * (ffct(n1, 5) / 2) * (mu_H[19] - mu_F[1] * mu_F[3]) + # Banner
      1 * (ffct(n1, 4) / 4) * (mu_H[2] - mu_F[1] * mu_F[3]) + # Diamond
      4 * (ffct(n1, 4) / 8) * (mu_H[8] - mu_F[1] * mu_F[3])) / ((ffct(n1, 4) / 8) * choose(n1, 2))) # Square
  sgm1[2, 3] <- sgm1[3, 2] <- sum((
    1 * (ffct(n1, 6) / 4) * (mu_H[6] - mu_F[2] * mu_F[3]) + # tsg2
      3 * (ffct(n1, 5) / 12) * (mu_H[5] - mu_F[2] * mu_F[3]) + # tsg3
      1 * (ffct(n1, 5) / 2) * (mu_H[4] - mu_F[2] * mu_F[3]) + # tsg4
      2 * (ffct(n1, 4) / 4) * (mu_H[2] - mu_F[2] * mu_F[3])) / ((ffct(n1, 4) / 8) * choose(n1, 3))) # tsg5
  #
  sgm2[1, 1] <- sum((2 * (ffct(n2, 3) / 2) * (mu_H[16] - mu_F[1]^2) + # Chery
    1 * choose(n2, 2) * (mu_H[17] - mu_F[1]^2)) / (choose(n2, 2)^2)) # Edge
  sgm2[2, 2] <- sum((2 * (ffct(n2, 5) / 8) * (mu_H[3] - mu_F[2]^2) + # Butterfly; tg1
    2 * (ffct(n2, 4) / 4) * (mu_H[2] - mu_F[2]^2) + # Diamond; tg2
    choose(n2, 3) * (mu_H[1] - mu_F[2]^2)) / (choose(n2, 3)^2)) # Triangle; tg0
  sgm2[3, 3] <- sum((2 * (ffct(n2, 7) / 8) * (mu_H[15] - mu_F[3]^2) + # sg1 (see end comments)
    6 * (ffct(n2, 6) / 48) * (mu_H[12] - mu_F[3]^2) + # sg2
    2 * (ffct(n2, 6) / 4) * (mu_H[13] - mu_F[3]^2) + # sg3
    2 * (ffct(n2, 6) / 4) * (mu_H[11] - mu_F[3]^2) + # sg4
    2 * (ffct(n2, 5) / 2) * (mu_H[9] - mu_F[3]^2) + # sg5
    6 * (ffct(n2, 5) / 12) * (mu_H[10] - mu_F[3]^2) + # sg6
    6 * (ffct(n2, 4) / 24) * (mu_H[7] - mu_F[3]^2) + # sg7
    1 * (ffct(n2, 4) / 8) * (mu_H[8] - mu_F[3]^2)) / (ffct(n2, 4) / 8)^2) #  Square
  sgm2[1, 2] <- sgm2[2, 1] <- sum((
    1 * (ffct(n2, 4) / 2) * (mu_H[18] - mu_F[1] * mu_F[2]) + # Shovel
      3 * choose(n2, 3) * (mu_H[1] - mu_F[1] * mu_F[2])) / (choose(n2, 3) * choose(n2, 2))) # Triangle
  sgm2[1, 3] <- sgm2[3, 1] <- sum((
    1 * (ffct(n2, 5) / 2) * (mu_H[19] - mu_F[1] * mu_F[3]) + # Banner
      1 * (ffct(n2, 4) / 4) * (mu_H[2] - mu_F[1] * mu_F[3]) + # Diamond
      4 * (ffct(n2, 4) / 8) * (mu_H[8] - mu_F[1] * mu_F[3])) / ((ffct(n2, 4) / 8) * choose(n2, 2))) # Square
  sgm2[2, 3] <- sgm2[3, 2] <- sum((
    1 * (ffct(n2, 6) / 4) * (mu_H[6] - mu_F[2] * mu_F[3]) + # tsg2
      3 * (ffct(n2, 5) / 12) * (mu_H[5] - mu_F[2] * mu_F[3]) + # tsg3
      1 * (ffct(n2, 5) / 2) * (mu_H[4] - mu_F[2] * mu_F[3]) + # tsg4
      2 * (ffct(n2, 4) / 4) * (mu_H[2] - mu_F[2] * mu_F[3])) / ((ffct(n2, 4) / 8) * choose(n2, 3))) # tsg5

  # Computing the p-value
  ######################
  g1 <- length(G1)
  g2 <- length(G2)
  hat_mu <- sqrt((g1 * g2) / (g1 + g2)) * (colMeans(mu_F1) - colMeans(mu_F2))
  hat_sgm <- g2 / (g1 + g2) * sgm1 / g1 + g1 / (g1 + g2) * sgm2 / g2
  statistic <- hat_mu %*% solve(hat_sgm) %*% hat_mu
  p.value <- 1 - pchisq(statistic, 3)

  # Output
  #######
  return(list(
    hat_mu = hat_mu,
    hat_sgm = hat_sgm,
    p.value = p.value,
    p.value_S = pnorm(hat_mu[3], 0, sqrt(hat_sgm[3, 3])),
    mu_F1 = mu_F1,
    mu_F2 = mu_F2,
    statistic = statistic
  ))
}
