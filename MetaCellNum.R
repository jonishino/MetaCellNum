library(doSNOW)
library(foreach)

#----------- Log-likelihood for one SNP data (D1, m1, D2, m2)
LL_eachSNP <- function(Nb,D1,m1,D2,m2,min_m1=2,p_bin,f_p1,pur1=1,pur2=1){
  if (length(f_p1)==1) f_p1 <- numeric(length(p_bin))+1/length(p_bin)
  p_plus_Mb <- cbind(p_bin*2, matrix(rep(0:Nb,length(p_bin)),nrow=length(p_bin),byrow=T))
  Mb_pr <- t(apply(p_plus_Mb,1, function(x) dbinom(x[-1], Nb, x[1])))
  pr <- ((f_p1 * dbinom(m1, D1, p_bin*pur1)) %*% (Mb_pr %*% dbinom(m2, D2, pur2*0:Nb/(2*Nb))))
  pr_min_m1 <- sum(apply(matrix(min_m1:D1), 1, function(m1_mat) (f_p1 * dbinom(m1_mat[1], D1, pur1*p_bin))))
  return(log(pr/pr_min_m1))
}

#----------- Log-likelihood for all SNP data (D1, m1, D2, m2)
LL_allSNPs <- function(Nb, dsin,min_m1,f_p1,p_bin,pur1,pur2){
   sum(apply(as.matrix(dsin), 1, function(x) LL_eachSNP(Nb=Nb,D1=x[1],m1=x[2],D2=x[3],m2=x[4],min_m1=min_m1,f_p1=f_p1,p_bin,pur1=pur1,pur2=pur2)))
}

#----------- Log-likelihood for each SNP data (D1, m1, D2, m2)
Disp_LL_allSNPs <- function(Nb, dsin,min_m1,f_p1,p_bin,pur1,pur2){
   apply(as.matrix(dsin), 1, function(x) LL_eachSNP(Nb=Nb,D1=x[1],m1=x[2],D2=x[3],m2=x[4],min_m1=min_m1,f_p1=f_p1,p_bin,pur1=pur1,pur2=pur2))
}

#----------- Search MLE of Nb for dataset of (D1, m1, D2, m2) by golden section search
mleNb <- function(dsin,min_m1=2,pur1,pur2,w=0.1,p_binN=100,p_binMin=10^-6,p_bin=NA,f_p1=NA,low_bd=1,upr_bd=1000,all_search_range=1:10, tolerance=3, disp=F, maxLLs=F){
    if(is.na(p_bin)[1] & is.na(f_p1)[1]){
        int_log = (log(0.5)-log(p_binMin))/(p_binN-1) #log-scale interval between i-th and (i+1)-th bin
        p_bin = exp(log(p_binMin) + int_log*0:(p_binN-1))
        cum_p1 = 1/p_bin - 1/0.5 #Williams et al.(2016) formulae (7)
        f_p1 = cum_p1[1:(p_binN-1)] - cum_p1[2:p_binN] 
        f_p1[p_binN] = sum(f_p1)/(1-w)*w
        f_p1 <- f_p1 /sum(f_p1)
    }
    else if(length(p_bin)==length(f_p1)){
        f_p1[p_binN] = sum(f_p1)/(1-w)*w
        f_p1 <- f_p1 /sum(f_p1)
    }
    else cat("Error in setting for p_bin and/or f_p1\n")
    gold <- 2/(sqrt(5)+1)
    x1 <- upr_bd-round(gold*(upr_bd-low_bd))
    x2 <- low_bd+round(gold*(upr_bd-low_bd))
        f1 <- LL_allSNPs(dsin=dsin,min_m1=min_m1,Nb=x1,f_p1=f_p1,p_bin=p_bin,pur1=pur1,pur2=pur2)
        f2 <- LL_allSNPs(dsin=dsin,min_m1=min_m1,Nb=x2,f_p1=f_p1,p_bin=p_bin,pur1=pur1,pur2=pur2)
    itr = 0
    while (abs(upr_bd - low_bd) > tolerance){
      itr <- itr+1
      if(disp==T) cat(itr,': Lower/Upper Bounds =', low_bd, '/', upr_bd, '\n')
      if (f2 < f1){
         upr_bd = x2
         x2 = x1
         f2 = f1
         x1 = upr_bd - round(gold*(upr_bd - low_bd))
         f1 = LL_allSNPs(dsin=dsin,min_m1=min_m1,Nb=x1,f_p1=f_p1,p_bin=p_bin,pur1=pur1,pur2=pur2)
      }
      else{
         low_bd = x1
         x1 = x2
         f1 = f2
         x2 = low_bd + round(gold*(upr_bd - low_bd))
         f2 = LL_allSNPs(dsin=dsin,min_m1=min_m1,Nb=x2,f_p1=f_p1,p_bin=p_bin,pur1=pur1,pur2=pur2)
      }
   }
   if(disp==T) cat('Final Lower/Upper Bounds =', low_bd, '/', upr_bd, '\n')
   if(disp==T) cat('\n......\n')
   bound.range <- low_bd:upr_bd
   LL<-c()
   for(i in bound.range) LL <- c(LL, LL_allSNPs(dsin=dsin,min_m1=min_m1,Nb=i,f_p1=f_p1,p_bin=p_bin,pur1=pur1,pur2=pur2))
   if(disp==T) cat("Log Likelihoods for Nb within Lower/Upper Bounds", LL, '\n')
   mle <- bound.range[which.max(LL)]
   LL_at_mle <- max(LL)
   if(disp==T) cat('MLE of Nb = ', mle,"(candidate)\n")
   if (length(all_search_range)>=2){
      LLs_for_all_search_range <- numeric(length(all_search_range))
      if(disp==T) cat('\n......\n')
      for (i in all_search_range) LLs_for_all_search_range[i] <- LL_allSNPs(dsin=dsin,min_m1=min_m1,Nb=i,f_p1=f_p1,p_bin=p_bin,pur1=pur1,pur2=pur2)
      LL_at_mle2 <- max(LLs_for_all_search_range)
      mle2 <- all_search_range[which.max(LLs_for_all_search_range)]
      if (LL_at_mle < LL_at_mle2) {
          if(disp==T) cat('Replacing MLE=', mle, "with",mle2, "due to complete search for Nb =", all_search_range, "with\nLog-Likelihoods of", LLs_for_all_search_range, '\n')
          mle <- mle2
      }
   }
   if(disp==T) cat('\n\nMLE of Nb = ', mle,"(final)\n\n\n")
   names(mle) <- "MLE"
   if(maxLLs==T){
#	mll <- max(c(LL_at_mle,LL_at_mle2))
	mll <- Disp_LL_allSNPs(dsin=dsin,min_m1=min_m1,Nb=x2,f_p1=f_p1,p_bin=p_bin,pur1=pur1,pur2=pur2)
#	names(mll) <- "MLL"
	return(list(MLE=mle, LLs=mll))
   }
   else return(mle)
}


#————— Nonparametric bootstrapping
mleNb_nbs <- function(dsin,pur1=1,pur2=1,w=0.1,p_binN=100,p_binMin=10^-6,p_bin=NA,f_p1=NA,low_bd=1,upr_bd=1000,min_m1=2,all_search_range=1:10, tolerance=3, disp=F,reps=100,clusters=8){
    registerDoSNOW(makeCluster(clusters, type="SOCK"))
    aaa <- foreach(i=1:reps, .combine = "c",.export=c("mleNb","LL_allSNPs","LL_eachSNP")) %dopar% {
       mleNb(dsin=dsin[sample(1:length(dsin$D1),replace=T),],pur1=pur1,pur2=pur2,w=w,min_m1=min_m1,p_binN=p_binN,p_binMin=p_binMin,p_bin=p_bin,f_p1=f_p1,low_bd=low_bd,upr_bd=upr_bd,all_search_range=all_search_range, tolerance=tolerance, disp=disp)
   }
   return(aaa)
}


