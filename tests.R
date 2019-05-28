require(igraph)
require(mixtools)
require(SDMTools)
require(igraph)
require(Matrix)
require(grid)
require(gridExtra)
library(ggplot2)
source("utilities.R")
source("network_sample_tests.R")
source("subgraph_counting.R")
pie <- read.csv("http://oeis.org/A000796/b000796.txt", header = FALSE, sep = " ")
set.seed(1)

##########
# SBM_CR
##########
rG <- function(n, B, P) { # Sample a graph uniformly over FSBM class with random size
  x <- runif(n) # Latent x_i s
  return(get.adjacency(sample_sbm(n, B, c(sum(x < P[1]), sum(x > P[1]))), sparse = F))
}

B1 <- matrix(c(.06, .02, .02, .06), 2, 2)
B2 <- matrix(c(.06, .02, .02, .05), 2, 2)
P <- c(.25, .75)
G <- sapply(pie[1:100, 2] + 30, rG, B = B1, P = P, simplify = F)
F1 <- SBm_CR(G, B1, P)
F2 <- SBm_CR(G, B2, P)
##########
# FSBm_CR0
##########
Bb <- function(k, u, v) {
  B <- matrix(v, k, k)
  diag(B) <- u
  B
}
rG <- function(n, B) {
  # Sample a graph uniformly over FSBM class with random size
  a <- c(0, sample_dirichlet(1, c(10, rep(1, nrow(B) - 1)))) # Probability of being in each group
  a <- cumsum(a) # Intervals
  x <- runif(n) # Latent x_i s
  P <- array(NA, nrow(B)) # Block sizes
  for (i in 1:nrow(B)) {
    P[i] <- sum((a[i] < x) * (x < a[i + 1]))
  }
  return(get.adjacency(sample_sbm(n, B, P), sparse = F))
}
B1 <- Bb(3, .06, .02)
B2 <- Bb(3, .055, .04)
G <- sapply(pie[1:200, 2] + 30, rG, B = B1, simplify = F)
F1 <- FSBm_CR0(G, B1)
F2 <- FSBm_CR0(G, B2)
##########
# FSBm_CR
##########
Bb <- function(k, u, v) {
  B <- matrix(v, k, k)
  diag(B) <- u
  return(B)
}
rG <- function(n, B) {
  # Sample a graph uniformly over FSBM class with random size
  a <- c(0, sample_dirichlet(1, c(10, rep(1, nrow(B) - 1)))) # Probability of being in each group
  a <- cumsum(a) # Intervals
  x <- runif(n) # Latent x_i s
  P <- array(NA, nrow(B)) # Block sizes
  for (i in 1:nrow(B)) {
    P[i] <- sum((a[i] < x) * (x < a[i + 1]))
  }
  return(get.adjacency(sample_sbm(n, B, P), sparse = F))
}

set.seed(1)
B1 <- Bb(3, .7, .2)
B2 <- Bb(5, .8, .55)
G <- sapply(pie[1:1000, 2] + 50, rG, B = B1, simplify = F)
F1 <- FSBm_CR(G, B1, do.plot = TRUE)
F2 <- FSBm_CR(G, B2, do.plot = TRUE)
############
# TSMPL_CR
############
rG <- function(n, B, P) {
  # Sample a graph uniformly over FSBM class with random size
  x <- runif(n) # Latent x_i s
  return(get.adjacency(sample_sbm(
    n,
    B,
    c(sum(x < P[1]), sum(x > P[1]))
  ),
  sparse = F
  ))
}
s_g <- function(A, p) {
  # Keep a fraction p of edges in A, unifornly at random
  as_adj(sample_gnp(nrow(A), p), sparse = F) * A
}
randomized_test <- function(G, p = 1) {
  G1 <- list()
  G2 <- list()
  for (i in 1:length(G)) {
    if (runif(1) > .5) {
      G1[[length(G1) + 1]] <- s_g(G[[i]], p)
    }
    else {
      G2[[length(G2) + 1]] <- s_g(G[[i]], p)
    }
  }
  TSMPL_CR(G1, G2)$p.value_S
}

# Test1
B1 <- matrix(c(.6, .2, .2, .1), 2, 2)
P1 <- c(.25, .75)
G1 <- sapply(pie[1:20, 2] + 20, rG, B = B1, P = P1, simplify = F)
B2 <- matrix(c(.9, .2, .2, .1), 2, 2)
P2 <- c(.25, .75)
G2 <- sapply(pie[101:128, 2] + 22, rG, B = B2, P = P2, simplify = F)
TSMPL_CR(G1, G2)

# Test2
B <- matrix(c(.5, .1, .1, .6), 2, 2)
P <- c(.5, .5)
G <- sapply(pie[1:200, 2] + 60, rG, B = B, P = P, simplify = F)
p.vals <- replicate(100, randomized_test(G, .05))
ks.test(p.vals, punif)
