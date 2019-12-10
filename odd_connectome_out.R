# Loading in the requiered sources
require(igraph)
source("network_sample_tests.R")
source("utilities.R")
source("subgraph_counting.R")

# Loading in the graphs
load(url("http://www.cis.jhu.edu/~parky/TT/Data/TT-DSDS01216-glist114-raw-LCCFALSE.rda"))

# Extracting the outlier
outlier_graph <- list(as_adjacency_matrix(glist[[107]]))
graph_sample <- lapply(glist[c(1:49,51:57)], as_adjacency_matrix)
rm(list = "glist")

# Computing the counts with TSMPL_CR_novar
test_output = TSMPL_CR_novar(graph_sample, outlier_graph)

# Conducting the test
set.seed(1)
n = nrow(test_output$mu_F1)
pvals = c() # pvalue while including the outlier
pvalsa = c() # pvalue while excluding the outlier
S = cov(test_output$mu_F1)
N = 1000
for (i in 1:N){
  sm = sample(n,n/2)
  m1 = colMeans(test_output$mu_F1[sm,])
  m2 = colMeans(rbind(test_output$mu_F1,test_output$mu_F2)[-(sm),])
  m2a = colMeans(test_output$mu_F1[-(sm),])
  pvals = c(pvals, 1- pchisq(((n+1)/2) * (m2-m1) %*% solve(2 * S) %*% (m2-m1), 3))
  pvalsa = c(pvalsa, 1- pchisq(((n+1)/2) * (m2a-m1) %*% solve(2 * S) %*% (m2a-m1), 3))
}
pval = min(pvals) * N
pvala = min(pvalsa) * N

(pval <0.05) * (pvala > 0.05)
