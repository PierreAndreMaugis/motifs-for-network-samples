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

m <- 15
k_max <- 4
rho <- 1/20/k_max
rho_bg <- .00
n <- 4*factorial(k_max)
model_parameters <-list()
for (k in 1:k_max){
  model_parameters[[k]] <- list(B=diag(k) * rho * k + rho_bg,
                                P=rep(1,k)/k)
}

R <- 5000
Pvals <- array(rep(1,k_max^2*R), dim=c(k_max, k_max, R))
for (i in 1:k_max){
  for (j in 1:k_max){
    for (t in 1:R){
      network_sample <- sapply(rep(n,m),
                               function(n, B, P){
                                 return(get.adjacency(sample_sbm(n,B,n*P),sparse=F))
                               },
                               B=model_parameters[[i]]$B,
                               P=model_parameters[[i]]$P,
                               simplify=F)
      ijt_test <- SBm_CR(network_sample,
                         B=model_parameters[[j]]$B,
                         P=model_parameters[[j]]$P,
                         do.plot=FALSE)
      Pvals[i,j,t] <- ijt_test$p.value
      }
  }
}
apply(Pvals,1:2,function(x){sum(x<.05)/length(x)})[-1,-1]
