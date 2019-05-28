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

##########
# Example 1
##########
pie <- read.csv("http://oeis.org/A000796/b000796.txt", header = FALSE, sep = " ")
rG <- function(n, B, P) { # Sample a graph uniformly over FSBM class with random size
  x <- runif(n) # Latent x_i s
  return(get.adjacency(sample_sbm(n, B, c(sum(x < P[1]), sum(x > P[1]))), sparse = F))
}
#
set.seed(1)
B1 <- matrix(c(.06, .02, .02, .06), 2, 2)
B2 <- matrix(c(.06, .02, .02, .05), 2, 2)
P <- c(.25, .75)
G <- sapply(pie[1:100, 2] + 30, rG, B = B1, P = P, simplify = F)
F1 <- SBm_CR(G, B1, P, do.plot = FALSE)
F2 <- SBm_CR(G, B2, P, do.plot = FALSE)
#
gl <- list()
labs <- c("C2", "C3", "C4")
for (i in 3:1) {
  j <- (3:1)[i]
  v <- ggplot() + theme_gray() + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
  y_range <- pmax(range(
    F1$conf_region[[j]][, 2],
    F1$hat_mu_F[-i][2],
    F2$conf_region[[j]][, 2],
    F2$hat_mu_F[-i][2]
  ) * c(.95, 1.05), 0)
  v <- v + scale_y_continuous(
    breaks = seq(min(y_range),
      max(y_range),
      length.out = 5
    ),
    limits = y_range
  )
  x_range <- pmax(range(
    F1$conf_region[[j]][, 1],
    F1$hat_mu_F[-i][1],
    F2$conf_region[[j]][, 1],
    F2$hat_mu_F[-i][1]
  ) * c(.95, 1.05), 0)
  v <- v + scale_x_continuous(
    breaks = seq(min(x_range),
      max(x_range),
      length.out = 5
    ),
    limits = x_range
  )
  v <- v + geom_polygon(
    data = data.frame(
      x = F1$conf_region[[j]][, 1],
      y = F1$conf_region[[j]][, 2]
    ),
    aes(x, y), fill = "black", alpha = .9, colour = NA
  )
  v <- v + geom_polygon(
    data = data.frame(
      x = F2$conf_region[[j]][, 1],
      y = F2$conf_region[[j]][, 2]
    ),
    aes(x, y), fill = "black", alpha = .5, colour = NA
  )
  v <- v + geom_point(
    data = data.frame(
      x = F1$mu_F[-i][1],
      y = F1$mu_F[-i][2]
    ),
    aes(x, y), color = "white"
  )
  v <- v + geom_point(
    data = data.frame(
      x = F2$mu_F[-i][1],
      y = F2$mu_F[-i][2]
    ),
    aes(x, y), color = "white"
  )
  v <- v + geom_point(
    data = data.frame(
      x = F1$hat_mu[-i][1],
      y = F1$hat_mu[-i][2]
    ),
    aes(x, y), color = "white", shape = 3, size = 2
  )
  v <- v + geom_point(
    data = data.frame(
      x = F1$hat_mu[-i][1],
      y = F1$hat_mu[-i][2]
    ),
    aes(x, y), color = "white", shape = 1, size = 2
  )
  v <- v + xlab((labs[-i])[1])
  v <- v + ylab((labs[-i])[2])
  gl[[i]] <- v
}
ggsave(
  filename = "Example1.pdf",
  plot = grid.arrange(gl[[3]], gl[[2]], gl[[1]], nrow = 1, ncol = 3),
  width = 12, height = 4, units = "in"
)
##########
# Example 2
##########
Bb <- function(k, u, v) {
  B <- matrix(v, k, k)
  diag(B) <- u
  B
}
pie <- read.csv("http://oeis.org/A000796/b000796.txt", header = FALSE, sep = " ")
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
B1 <- Bb(3, .06, .02)
B2 <- Bb(3, .055, .04)
G <- sapply(pie[1:200, 2] + 30, rG, B = B1, simplify = F)
F1 <- FSBm_CR0(G, B1, do.plot = FALSE)
F2 <- FSBm_CR0(G, B2, do.plot = FALSE)
level <- .95
gl <- list()
labs <- c("C2", "C3", "C4")
for (i in 3:1) {
  j <- (3:1)[i]
  v <- ggplot() + theme_gray() + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
  y_range <- pmax(range(
    F1$conf_region[[j]][-1, 2],
    F1$hat_mu_F[-i][2],
    F2$conf_region[[j]][-1, 2],
    F2$hat_mu_F[-i][2]
  ) * c(.99, 1.01), 0)
  v <- v + scale_y_continuous(
    breaks = seq(min(y_range),
      max(y_range),
      length.out = 5
    ),
    limits = y_range
  )
  x_range <- pmax(range(
    F1$conf_region[[j]][-1, 1],
    F1$hat_mu_F[-i][1],
    F2$conf_region[[j]][-1, 1],
    F2$hat_mu_F[-i][1]
  ) * c(.99, 1.01), 0)
  v <- v + scale_x_continuous(
    breaks = seq(min(x_range),
      max(x_range),
      length.out = 5
    ),
    limits = x_range
  )
  super_ellipse_1 <- c()
  for (l in 1:length(F1$u)) {
    super_ellipse_1 <- rbind(
      super_ellipse_1,
      ellipse(F1$mu_F[l, -i],
        F1$sgm[-i, -i, l] / length(G)^2,
        alpha = 1 - level,
        draw = F
      )
    )
  }
  super_ellipse_contour_1 <- super_ellipse_1[chull(super_ellipse_1), ]
  v <- v + geom_polygon(
    data = data.frame(
      x = super_ellipse_contour_1[, 1],
      y = super_ellipse_contour_1[, 2]
    ),
    aes(x, y),
    fill = "gray60",
    alpha = .5,
    colour = "gray50",
    linetype = "dotted"
  )
  inner_ellipse_1 <- emb_contour(F1$mu_F[, -i])
  v <- v + geom_polygon(
    data = data.frame(
      x = inner_ellipse_1[, 1],
      y = inner_ellipse_1[, 2]
    ),
    aes(x, y),
    fill = "gray50",
    alpha = 1,
    colour = "gray50"
  )
  super_ellipse_2 <- c()
  for (l in 1:length(F2$u)) {
    super_ellipse_2 <- rbind(
      super_ellipse_2,
      ellipse(F2$mu_F[l, -i],
        F2$sgm[-i, -i, l] / length(G)^2,
        alpha = 1 - level,
        draw = FALSE
      )
    )
  }
  super_ellipse_contour_2 <- super_ellipse_2[chull(super_ellipse_2), ]
  v <- v + geom_polygon(
    data = data.frame(
      x = super_ellipse_contour_2[, 1],
      y = super_ellipse_contour_2[, 2]
    ),
    aes(x, y),
    fill = "black",
    alpha = .25,
    colour = "gray20",
    linetype = "dashed"
  )
  inner_ellipse_2 <- emb_contour(F2$mu_F[, -i])
  v <- v + geom_polygon(
    data = data.frame(
      x = inner_ellipse_2[, 1],
      y = inner_ellipse_2[, 2]
    ),
    aes(x, y),
    fill = "black",
    alpha = 1,
    colour = "black"
  )
  v <- v + geom_point(
    data = data.frame(
      x = F1$hat_mu[-i][1],
      y = F1$hat_mu[-i][2]
    ),
    aes(x, y), color = "white", shape = 3, size = 2
  )
  v <- v + geom_point(
    data = data.frame(
      x = F1$hat_mu[-i][1],
      y = F1$hat_mu[-i][2]
    ),
    aes(x, y), color = "white", shape = 1, size = 2
  )
  v <- v + xlab((labs[-i])[1])
  v <- v + ylab((labs[-i])[2])
  gl[[i]] <- v
}
ggsave(
  filename = "Example2.pdf",
  plot = grid.arrange(gl[[3]], gl[[2]], gl[[1]], nrow = 1, ncol = 3),
  width = 12, height = 4, units = "in"
)
##########
# Example 3
##########
Bb <- function(k, u, v) {
  B <- matrix(v, k, k)
  diag(B) <- u
  return(B)
}
pie <- read.csv("http://oeis.org/A000796/b000796.txt", header = FALSE, sep = " ")
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
gl <- list()
labs <- c("C2", "C3", "C4")
for (i in 3:1) {
  j <- (3:1)[i]
  v <- ggplot() + theme_gray() + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
  y_range <- pmax(range(
    F1$aux_sph[[j]][, 2],
    F1$hat_mu_F[-i][2],
    F2$aux_sph[[j]][, 2],
    F2$hat_mu_F[-i][2]
  ) * c(.99, 1.01), 0)
  v <- v + scale_y_continuous(
    breaks = seq(min(y_range),
      max(y_range),
      length.out = 5
    ),
    limits = y_range
  )
  x_range <- pmax(range(
    F1$aux_sph[[j]][, 1],
    F1$hat_mu_F[-i][1],
    F2$aux_sph[[j]][, 1],
    F2$hat_mu_F[-i][1]
  ) * c(.99, 1.01), 0)
  v <- v + scale_x_continuous(
    breaks = seq(min(x_range),
      max(x_range),
      length.out = 5
    ),
    limits = x_range
  )
  v <- v + geom_polygon(
    data = data.frame(
      x = F1$aux_sph[[j]][, 1],
      y = F1$aux_sph[[j]][, 2]
    ),
    aes(x, y),
    fill = "gray60",
    alpha = .5,
    colour = "gray50",
    linetype = "dotted"
  )
  v <- v + geom_polygon(
    data = data.frame(
      x = F2$aux_sph[[j]][, 1],
      y = F2$aux_sph[[j]][, 2]
    ),
    aes(x, y),
    fill = "black",
    alpha = .5,
    colour = "gray20",
    linetype = "dashed"
  )
  v <- v + geom_polygon(
    data = data.frame(
      x = F1$ch0[[j]][, 1],
      y = F1$ch0[[j]][, 2]
    ),
    aes(x, y),
    fill = "gray60",
    alpha = .5,
    colour = "gray50"
  )
  v <- v + geom_polygon(
    data = data.frame(
      x = F2$ch0[[j]][, 1],
      y = F2$ch0[[j]][, 2]
    ),
    aes(x, y),
    fill = "black",
    alpha = .5,
    colour = "gray20"
  )
  # polygon(F1$ch0[[j]],border=1,lwd=2)
  # polygon(F2$ch0[[j]],border=1,lwd=2)
  v <- v + geom_point(
    data = data.frame(
      x = F1$hat_mu[-i][1],
      y = F1$hat_mu[-i][2]
    ),
    aes(x, y), color = "white", shape = 3, size = 2
  )
  v <- v + geom_point(
    data = data.frame(
      x = F1$hat_mu[-i][1],
      y = F1$hat_mu[-i][2]
    ),
    aes(x, y), color = "white", shape = 1, size = 2
  )
  v <- v + xlab((labs[-i])[1])
  v <- v + ylab((labs[-i])[2])
  gl[[i]] <- v
}
ggsave(
  filename = "Example3.pdf",
  plot = grid.arrange(gl[[3]], gl[[2]], gl[[1]], nrow = 1, ncol = 3),
  width = 12, height = 4, units = "in"
)

############
# Finite sample experiments
############
# Set-Up
rB <- function(k) {
  B <- matrix(runif(k^2), k, k)
  B[lower.tri(B)] <- B[upper.tri(B)]
  B
}
rP <- function(k) {
  sample_dirichlet(1, rep(1, k))
}
rG <- function(n, B, P) { # Sample a graph
  x <- runif(n) # Latent x_i s
  a <- array(NA, nrow(B)) # Block sizes
  P2 <- c(0, cumsum(P))
  for (i in 1:nrow(B)) {
    a[i] <- sum((P2[i] < x) * (x < P2[i + 1]))
  }
  get.adjacency(sample_sbm(n, B, a))
}
rTst <- function(n, k) {
  B <- rB(k)
  P <- rP(k)
  G <- lapply(n, rG, B = B, P = P)
  Tst <- SBm_CR(G, B, P, do.plot = F)
  return(c(Tst$p.value, Tst$dst))
}
# Simulation Loop
R <- 500
N <- c(100, 200, 400, 800, 1600, 3200)
pie <- read.csv("http://oeis.org/A000796/b000796.txt",
  header = FALSE, sep = " "
)
n <- pie[1:max(N), 2] + 8
Tsts <- array(NA, c(2, R, length(N)))
ptm <- proc.time()
for (i in 1:length(N)) {
  Tsts[, , i] <- replicate(R, rTst(n[1:N[i]], 2))
  save(list = "Tsts", file = ("~/Desktop/Tsts.Rmat"))
  cat(paste0(
    round(100 * i / length(N), 2), "% ",
    round((proc.time() - ptm)[["elapsed"]] / 60, 2),
    " mns elapsed at ", format(Sys.time(), "%X"), "
"
  ))
  flush.console()
}

# Plotting
quartz()
par(mfrow = c(1, 2))
df1 <- data.frame(y = c(t(Tsts[2, , ])), x = factor(rep(1:length(N))))
boxplot(log(y) ~ x, data = df1, axes = F, xlab = "Network Sample Size (N)", ylab = "Root Squares Error (log scale)")
y_axis <- seq(par("usr")[3], par("usr")[4], length.out = 5)
regline <- lm(log(y) ~ as.numeric(x), data = df1)
abline(regline, col = 2, lwd = 2, lty = "solid")
options("scipen" = -100, "digits" = 2)
axis(2, at = y_axis, labels = paste(signif(exp(y_axis), 2)))
options("scipen" = 100, "digits" = 2)
axis(1, at = 1:length(N), labels = N)

abline(h = .1, col = 2, lwd = 2, lty = "solid")
qqplot(qunif(seq(0, 1, by = .01)), Tsts[1, , ],
  axes = F,
  xlab = "Uniform Distribution",
  ylab = "p-values"
)
qqline(Tsts[1, , ], distribution = qunif)
axis(1, at = seq(0, 1, by = .2))
axis(2, at = seq(0, 1, by = .2))
