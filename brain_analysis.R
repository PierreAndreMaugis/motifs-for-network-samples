############
# Brain network data analysis
############
con <- url("http://www.cis.jhu.edu/~parky/MRN/fibergraph.Rbin")
load(con)
close(con)
con <- url("http://www.cis.jhu.edu/~parky/MRN/newcci.txt")
cci <- scan(con)
close(con)
G <- list()
for (i in 1:length(fibergraph.list)) {
  gi <- fibergraph.list[[i]] # Upper-triangle only
  n <- gi@Dim[1]
  gi <- (gi + t(gi)) # Aymmetrization, automatic hollow
  gi <- ifelse(gi > 0, 1, 0) # Binarization
  gi <- matrix(gi, nrow = n)
  G[[i]] <- gi
}
## Testing uniform subsampling
set.seed(2)
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
  return(TSMPL_CR(G1, G2)$p.value_S)
}
p.vals <- replicate(100, randomized_test(G))
ks.test(p.vals, punif) # D = 0.094574, p-value = 0.3327

## Analysis
covariates_test <- function(G, X, q = .5, p = 1) {
  G1 <- G2 <- list()
  q <- min(q, 1 - q)
  p <- min(1, p)
  for (i in 1:length(G)) {
    if (!is.na(X[i])) {
      if (X[i] <= quantile(X, q, na.rm = T)) {
        G1[[length(G1) + 1]] <- s_g(G[[i]], p)
      }
      else {
        if (X[i] > quantile(X, 1 - q, na.rm = T)) {
          G2[[length(G2) + 1]] <- s_g(G[[i]], p)
        }
      }
    }
  }
  TSMPL_CR(G1, G2)
}
ttt1 <- covariates_test(G, cci, .5)
ttt2 <- covariates_test(G, cci, .4)
ttt3 <- covariates_test(G, cci, .3)
ttt4 <- covariates_test(G, cci, .2)
ttt5 <- covariates_test(G, cci, .1) # Sample too small
#
Tbl <- 2 * rbind(
  pnorm(ttt1$hat_mu / sqrt(diag(ttt1$hat_sgm))),
  pnorm(ttt2$hat_mu / sqrt(diag(ttt2$hat_sgm))),
  pnorm(ttt3$hat_mu / sqrt(diag(ttt3$hat_sgm))),
  pnorm(ttt4$hat_mu / sqrt(diag(ttt4$hat_sgm))),
  pnorm(ttt5$hat_mu / sqrt(diag(ttt5$hat_sgm)))
)
colnames(Tbl) <- c("Edge", "Triangle", "Square")
rownames(Tbl) <- c(.5, .4, .3, .2, .1)
Tbl <- round(Tbl, 3)

cG <- lapply(G, function(x) {
  1 - x
})
tta1 <- covariates_test(cG, cci, .5)
tta2 <- covariates_test(cG, cci, .4)
tta3 <- covariates_test(cG, cci, .3)
tta4 <- covariates_test(cG, cci, .2)
tta5 <- covariates_test(cG, cci, .1) # Sample too small
#
Tbla <- 2 * (1 - rbind(
  pnorm(tta1$hat_mu / sqrt(diag(tta1$hat_sgm))),
  pnorm(tta2$hat_mu / sqrt(diag(tta2$hat_sgm))),
  pnorm(tta3$hat_mu / sqrt(diag(tta3$hat_sgm))),
  pnorm(tta4$hat_mu / sqrt(diag(tta4$hat_sgm))),
  pnorm(tta5$hat_mu / sqrt(diag(tta5$hat_sgm)))
))
colnames(Tbla) <- c("Edge", "Triangle", "Square")
rownames(Tbla) <- c(.5, .4, .3, .2, .1)
Tbla <- round(Tbla, 3)

#####################
# Replication of analysis using that graphs are iid and vertex matched
#
# Tbl2 <- matrix(c(
#  t.test(ttt1$mu_F1[,1],ttt1$mu_F2[,1],var.equal=T)$p.value,
#  t.test(ttt1$mu_F1[,2],ttt1$mu_F2[,2],var.equal=T)$p.value,
#  t.test(ttt1$mu_F1[,3],ttt1$mu_F2[,3],var.equal=T)$p.value,
#  t.test(ttt2$mu_F1[,1],ttt2$mu_F2[,1],var.equal=T)$p.value,
#  t.test(ttt2$mu_F1[,2],ttt2$mu_F2[,2],var.equal=T)$p.value,
#  t.test(ttt2$mu_F1[,3],ttt2$mu_F2[,3],var.equal=T)$p.value,
#  t.test(ttt3$mu_F1[,1],ttt3$mu_F2[,1],var.equal=T)$p.value,
#  t.test(ttt3$mu_F1[,2],ttt3$mu_F2[,2],var.equal=T)$p.value,
#  t.test(ttt3$mu_F1[,3],ttt3$mu_F2[,3],var.equal=T)$p.value,
#  t.test(ttt4$mu_F1[,1],ttt4$mu_F2[,1],var.equal=T)$p.value,
#  t.test(ttt4$mu_F1[,2],ttt4$mu_F2[,2],var.equal=T)$p.value,
#  t.test(ttt4$mu_F1[,3],ttt4$mu_F2[,3],var.equal=T)$p.value,
#  t.test(ttt5$mu_F1[,1],ttt5$mu_F2[,1],var.equal=T)$p.value,
#  t.test(ttt5$mu_F1[,2],ttt5$mu_F2[,2],var.equal=T)$p.value,
#  t.test(ttt5$mu_F1[,3],ttt5$mu_F2[,3],var.equal=T)$p.value),
#  ncol=3,byrow=T)
# colnames(Tbl2) <- c('Edge','Triangle','Square')
# rownames(Tbl2) <- c(.5,.4,.3,.2,.1)
# Tbl2 <- round(Tbl2,3);
#
# Tbla2 <- matrix(c(
#  t.test(tta1$mu_F1[,1],tta1$mu_F2[,1],var.equal=T)$p.value,
#  t.test(tta1$mu_F1[,2],tta1$mu_F2[,2],var.equal=T)$p.value,
#  t.test(tta1$mu_F1[,3],tta1$mu_F2[,3],var.equal=T)$p.value,
#  t.test(tta2$mu_F1[,1],tta2$mu_F2[,1],var.equal=T)$p.value,
#  t.test(tta2$mu_F1[,2],tta2$mu_F2[,2],var.equal=T)$p.value,
#  t.test(tta2$mu_F1[,3],tta2$mu_F2[,3],var.equal=T)$p.value,
#  t.test(tta3$mu_F1[,1],tta3$mu_F2[,1],var.equal=T)$p.value,
#  t.test(tta3$mu_F1[,2],tta3$mu_F2[,2],var.equal=T)$p.value,
#  t.test(tta3$mu_F1[,3],tta3$mu_F2[,3],var.equal=T)$p.value,
#  t.test(tta4$mu_F1[,1],tta4$mu_F2[,1],var.equal=T)$p.value,
#  t.test(tta4$mu_F1[,2],tta4$mu_F2[,2],var.equal=T)$p.value,
#  t.test(tta4$mu_F1[,3],tta4$mu_F2[,3],var.equal=T)$p.value,
#  t.test(tta5$mu_F1[,1],tta5$mu_F2[,1],var.equal=T)$p.value,
#  t.test(tta5$mu_F1[,2],tta5$mu_F2[,2],var.equal=T)$p.value,
#  t.test(tta5$mu_F1[,3],tta5$mu_F2[,3],var.equal=T)$p.value),
#  ncol=3,byrow=T)
# colnames(Tbla2) <- c('Edge','Triangle','Square');
# rownames(Tbla2) <- c(.5,.4,.3,.2,.1);
# Tbla2 <- round(Tbla2,3);
#
# Tbl;Tbl2;Tbla;Tbla2
# covariates_split <-function(X,q=.5){
#  X1 <- X2 <- list();
#  q <- min(q,1-q);
#  for (i in 1:length(X)) {
#    if (!is.na(X[i])) {
#      if (X[i]<=quantile(X,q,na.rm=T)) {
#        X1[[length(X1)+1]] <- X[i]}
#      else {if (X[i]>quantile(X,1-q,na.rm=T)) {
#        X2[[length(X2)+1]] <- X[i]}}
#    }}
#  return(c(unlist(X1),unlist(X2)))
# }
# cci2 <- covariates_split(cci,.5)
# mu_FF <- rbind(ttt1$mu_F1,ttt1$mu_F2)
# mu_FX <- rbind(tta1$mu_F1,tta1$mu_F2)
#
# cor(cbind(cci2,mu_FF,mu_FX))
#
# # Simple regression
# reg <- lm(cci2~mu_FF[,1]+mu_FF[,2]+mu_FF[,3]+
#             mu_FX[,2]+mu_FX[,3])
# summary(reg)
#
# # Quantile regression
# library('quantreg')
# x <- mu_FF[,2]
# qreg <- rq(cci2~mu_FX[,2]+x,tau=seq(.05,.95,by=.05))
# plot(summary(qreg), parm="x")
#
# plot(rq(cci2~mu_FF[,2]+mu_FX[,3],tau=seq(.05,.95,by=.05)))
#
# # Parametric method of moment test (second order null verification)
# source('utilities.R')
# MM_test <-function(G,X,q=.5,p=1){
#  G1 <- G2 <- list();
#  q <- min(q,1-q); p <- min(1,p);
#  for (i in 1:length(G)) {
#    if (!is.na(X[i])) {
#      if (X[i]<=quantile(X,q,na.rm=T)) {
#        G1[[length(G1)+1]] <- s_g(G[[i]],p)}
#      else {if (X[i]>quantile(X,1-q,na.rm=T)) {
#        G2[[length(G2)+1]] <- s_g(G[[i]],p)}}
#  }}
#  hat_B <- Bm_MM(G1)$hat_B
#  FSBm_CR(G2,hat_B)
# }
# X = MM_test(G,cci,.3)
#
