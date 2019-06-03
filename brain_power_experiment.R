# Largely based on http://www.cis.jhu.edu/~parky/MRN/omnicci.html
require(foreach)
con <- url("http://www.cis.jhu.edu/~parky/MRN/fibergraph.Rbin")
load(con)
close(con)
con <- url("http://www.cis.jhu.edu/~parky/MRN/newcci.txt")
cci <- scan(con)
close(con)
source("network_sample_tests.R")
source("utilities.R")
source("subgraph_counting.R")

glist <- list()
for (i in 1:length(fibergraph.list)) {
  gi <- fibergraph.list[[i]]; ## upper-triangle only
  n  <- gi@Dim[1];
  gi <- (gi + t(gi)); ## symmetrization, automatic hollow
  gi <- ifelse(gi>0,1,0); gi <- matrix(gi,nrow=n); ## binarization
  glist[[i]] <- graph_from_adjacency_matrix(gi, mode='undirected', diag=F);
}

empty <- which(sapply(glist, ecount)==0 | sapply(glist, vcount)!=70)
if (length(empty) > 0) {
  glist <- glist[-empty]
  cci <- cci[-empty]
} else {
  glist <- glist
}

which.na <- which(is.na(cci))
if (length(which.na) > 0) {
  cci <- cci[-which.na]
  glist <- glist[-which.na]
}
m <- length(cci)
glist <- lapply(glist, function(x) set_graph_attr(x, "cci", value=cci))

source(url("http://www.cis.jhu.edu/~parky/MRN/ccisim.R"))

thresh <- quantile(cci, c(.2,.8))
dhat <- 10
Svec <- case_when(
  cci <= thresh[1] ~ 1,
  cci >= thresh[2] ~ 2,
  TRUE ~ 0
)

gbinary <- glist
Alist <- lapply(gbinary, get.adjacency)
Abar <- Reduce('+', Alist) / length(Alist)
Xhat <- full.ase(Abar, dhat)$Xhat
Phat <- Xhat %*% t(Xhat)

Abar1 <- Reduce('+', Alist[Svec==1]) / sum(Svec==1)
Abar2 <- Reduce('+', Alist[Svec==2]) / sum(Svec==2)
Xhat1 <- full.ase(Abar1, dhat)$Xhat
Xhat2 <- full.ase(Abar2, dhat)$Xhat
Phat1 <- Xhat1 %*% t(Xhat1)
Phat2 <- Xhat2 %*% t(Xhat2)


runMC <- function(P0, P1, P2, ell, g0, gA, nmc)
{
  
  T0 <- foreach(mc=1:nmc,.combine='rbind') %dopar% {
    g01 <- replicate(ell,rg.sample(P0,P1,g0),simplify='list')
    g02 <- replicate(ell,rg.sample(P0,P2,g0),simplify='list')
    T0a <- calcT(Reduce('+',g01)/ell,Reduce('+',g02)/ell) # null (nonparT, semiparT)
    T0b <- TSMPL_CR(g01,g02)
    c(unlist(T0a),T0b$statistic)
  }
  rownames(T0) <- 1:nmc
  colnames(T0) <- c("nonpar","semipar","motif")#,"semipar.omni")
  
  TA <- NULL
  TA[[1]] <- cbind(1:mc, g0, T0)
  colnames(TA[[1]]) <- c("mc","gamma","nonpar","semipar",'motif')
  for (i in 2:length(gA)) {
    gamma <- gA[i]
    
    TA[[i]] <- foreach(mc=1:nmc,.combine='rbind') %dopar% {
      gA1 <- replicate(ell,rg.sample(P0,P1,gamma),simplify='list')
      gA2 <- replicate(ell,rg.sample(P0,P2,gamma),simplify='list')
      TAa <- calcT(Reduce('+',gA1)/ell,Reduce('+',gA2)/ell) # alt (nonparT, semiparT)
      TAb <- TSMPL_CR(gA1,gA2)
      c(mc=mc, gamma=gamma, unlist(TAa), TAb$statistic)
    }
    
    rownames(TA[[i]]) <- 1:nmc
    colnames(TA[[i]]) <- c("mc","gamma","nonpar","semipar",'motif')#,"semipar.omni")
  }
  
  return(list(T0=T0, TA=TA))
}

calcPwr <- function(T0, TA, alpha=0.05)
{
  nmc <- nrow(T0)
  cv.nonpar <- quantile(T0[,"nonpar"], probs = c(1-alpha))
  cv.semipar <- quantile(T0[,"semipar"], probs = c(1-alpha))
  cv.motif <- quantile(T0[,"motif"], probs = c(1-alpha))
  
  pwr.nonpar <- sum(TA[,"nonpar"] > cv.nonpar) / nmc
  pwr.semipar <- sum(TA[,"semipar"] > cv.semipar) / nmc
  pwr.motif <- sum(TA[,"motif"] > cv.motif) / nmc
  
  return(c(cv.nonpar=cv.nonpar, cv.semipar=cv.semipar, cv.motif=cv.motif,
           pwr.nonpar=pwr.nonpar, pwr.semipar=pwr.semipar, pwr.motif=pwr.motif))
}   

                
# Running the experiemnt (change the three parameters bellow to change how large you want the experiment to be) 
nmc <- 53
ell <- 23
g0 <- 0
gA <- seq(0, .5, by=.1ß)
Tlist <- runMC(Phat, Phat1, Phat2, ell, g0, gA, nmc)

# Plotting output
df.pwr <- as.data.frame(sapply(Tlist$TA, function(x) calcPwr(Tlist$T0, x)))
names(df.pwr) <- paste0("gamma = ", gA)

TA <- lapply(1:length(Tlist$TA), function(x) 
  cbind(Tlist$TA[[x]], 
        cv.nonpar=df.pwr["cv.nonpar.95%",x],
        cv.semipar=df.pwr["cv.semipar.95%",x],
        cv.motif=df.pwr["cv.motif.95%",x],
        pwr.nonpar=df.pwr["pwr.nonpar",x],
        pwr.semipar=df.pwr["pwr.semipar",x],
        pwr.motif=df.pwr["pwr.motif",x]))#,
TA <- as.data.frame(Reduce('rbind',TA))
names(TA)[3:5] <- paste0("stat.",names(TA)[3:5])
df.ta <- TA %>%
  gather(key, value, -mc, -gamma) %>%
  extract(key, c("stat", "test"), "(^.*)\\.(.*)") %>%
  spread(stat, value) 
df.ta$gamma <- factor(df.ta$gamma)

df.cv <- data.frame(test=c("nonpar","semipar","motif"), 
                    cv=c(df.pwr[1,1], df.pwr[2,1], df.pwr[3,1]))
df.pw <- data.frame(gamma=factor(gA), t(df.pwr[-c(1,2,3),]))
names(df.pw) <- c("gamma", "nonpar","semipar","motif")
df.pw <- melt(df.pw, variable.name = "test", value.name = "power")
df.pw <- df.pw %>% mutate(cv=rep(df.cv$cv, each=length(gA)))

p <- ggplot(df.ta, aes(x=stat, y=reorder(gamma, desc(gamma)), fill=gamma)) +
  facet_wrap(~test, ncol=3, scales="free_x") +
  stat_density_ridges(alpha=0.5, scale=1.5, adjust=2) + 
  theme(legend.position = "none") + labs(y="gamma") +
  geom_vline(data=df.cv, aes(xintercept=cv), linetype="dashed") 
p
p + geom_label(data = df.pw, aes(y=reorder(gamma,desc(gamma)), x=cv, label=sprintf("%1.2f",power)), size=3)

#ggplot(df.pw, aes(x=gamma, y=power, color=test, fill=test)) +
#  ylim(0,1) + 
#  geom_bar(stat="identity", position = "dodge", width=0.5, alpha=0.7)

pvales <- matrix(NA,length(gA)-1,2)
colnames(pvales) <- c('motif > semipar','nonpar > motif')
rownames(pvales) <- gA[-1]
for (i in 1:(length(gA)-1)){
  sset <- (TA["gamma"]==gA[i+1])
  pvales[i,1] <- t.test((TA[sset, "stat.motif"] > TA[sset,"cv.motif"])*1,
                        (TA[sset, "stat.semipar"] > TA[sset,"cv.semipar"])*1,
                        alternative="less", paiered=T)$p.value
  pvales[i,2] <- t.test((TA[sset, "stat.nonpar"] > TA[sset,"cv.nonpar"])*1,
                        (TA[sset, "stat.motif"] > TA[sset,"cv.motif"])*1,
                        alternative="less", paiered=T)$p.value
}
pvales>.05
