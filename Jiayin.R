Shihn <- function(mu1,mu2,V1,V2,n1,n2){
  wi <- (mu1 - mu2) / sqrt( V1 / n1 + V2 / n2 )
  mu_w <- mean(wi)
  nc <- ( (wi-mu_w) ^ 2 == max( (wi-mu_w) ^ 2 ) )
  out <- list(wi=wi,mu_w=mu_w,nc=nc)
  
  return(out)
}

ShihnSup <- function(mu1,mu2,V1,V2,n1,n2){
  wi <- (mu1 - mu2) / sqrt( V1 / n1 + V2 / n2 )
  mu_w <- mean(wi)
  nc <- ( (wi-mu_w) == min( (wi-mu_w) ) )
  out <- list(wi=wi,mu_w=mu_w,nc=nc)
  
  return(out)
}

Shihb <- function(p1,p2,n1,n2){
  # p <- (n1*p1+n2*p2) / (n1+n2)
  
  wi <- (p1-p2) / ( sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2) )
  mu_w <- mean(wi)
  nc <- ( (wi-mu_w) ^ 2 == max( (wi-mu_w) ^ 2 ) )
  out <- list(wi=wi,mu_w=mu_w,nc=nc)
  
  return(out)
}

ShihbSup <- function(p1,p2,n1,n2){
  # p <- (n1*p1+n2*p2) / (n1+n2)
  
  wi <- (p1-p2) / ( sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2) )
  mu_w <- mean(wi)
  nc <- ( (wi-mu_w) == min( (wi-mu_w) ) )
  out <- list(wi=wi,mu_w=mu_w,nc=nc)
  
  return(out)
}

var2g <- function(dataT,dataR){
  if (is.vector(dataT)){
    # C <- length(dataT)
    C1 <- sum(!is.na(dataT))
    C2 <- sum(!is.na(dataR))
    tempmat <- matrix(dataT[!is.na(dataT)],nrow=C1,ncol=C2,byrow=FALSE)-matrix(dataR[!is.na(dataR)],nrow=C1,ncol=C2,byrow=TRUE)
    var2g <- var(as.vector(tempmat))
  } else if (is.matrix(dataT)) { 
    L <- nrow(dataT)
    # C <- ncol(dataT)
    var2g <- rep(0,L)
    for (i in 1:L){
      C1 <- sum(!is.na(dataT[i,]))
      C2 <- sum(!is.na(dataR[i,]))
      tempmat <- matrix(dataT[i,!is.na(dataT[i,])],nrow=C1,ncol=C2,byrow=FALSE)-matrix(dataR[i,!is.na(dataR[i,])],nrow=C1,ncol=C2,byrow=TRUE)
      var2g[i] <- var(as.vector(tempmat))
    }
  } else {
    print("data typ error for var2g function!")
  }
  return(var2g)
}

Tse <- function(muS,muG,nS,nG,VS,VG,p0=0.8,d=0.2,b=0.1){
  K <- length(muS)
  muG <- rep(muG[1],K)
  nG <- rep(nG[1],K)
  VG <- rep(VG[1],K)
  
  hp <- pnorm( ( -log(1-d)-(muS-muG) )/sqrt(VS+VG) )-pnorm( ( log(1-d)-(muS-muG) )/sqrt(VS+VG) )
  z1 <- ( log(1-d)-(muS-muG) )/sqrt(VS+VG)
  z2 <- ( -log(1-d)-(muS-muG) )/sqrt(VS+VG)
  pmuS <- (-1/sqrt(VS+VG))*(dnorm(z2)-dnorm(z1))
  pmuG <- -pmuS
  pV <- (-1/2/(VS+VG))*(z2*dnorm(z2)-z1*dnorm(z1))
  pmu2 <- (-1/(VS+VG))*(z2*dnorm(z2)-z1*dnorm(z1))
  pV2 <- 1/4/((VS+VG)^2)*((3*z2-z2^3)*dnorm(z2)-(3*z1-z1^3)*dnorm(z1))
  Bp <- -pV*VS/nS-pV*VG/nG+1/2*(pmu2*VS/nS+pmu2*VG/nG+pV2*2*(VS^2)/nS+pV2*2*(VG^2)/nG)
  Cp <- pmuS^2*VS/nS+pmuG^2*VG/nG+pV^2*2*(VS^2)/nS+pV^2*2*(VG^2)/nG
  
  nc <- hp<=p0+Bp+qnorm(1-b)*sqrt(Cp)
  out <- list(hp=hp,Bp=Bp,Cp=Cp,nc=nc,d=d,p0=p0,b=b)
  
  return(out)
}

TseSup <- function(muS,muG,nS,nG,VS,VG,p0=0.8,d=0.2,b=0.1){
  K <- length(muS)
  muG <- rep(muG[1],K)
  nG <- rep(nG[1],K)
  VG <- rep(VG[1],K)
  
  hp <- 1-pnorm( ( log(1-d)-(muS-muG) )/sqrt(VS+VG) )
  z1 <- ( log(1-d)-(muS-muG) )/sqrt(VS+VG)
  pmuS <- 1/sqrt(VS+VG)*dnorm(z1)
  pmuG <- -pmuS
  pV <- 1/2/(VS+VG)*z1*dnorm(z1)
  pmu2 <- 1/(VS+VG)*z1*dnorm(z1)
  pV2 <- 1/4/((VS+VG)^2)*(-(3*z1-z1^3)*dnorm(z1))
  Bp <- -pV*VS/nS-pV*VG/nG+1/2*(pmu2*VS/nS+pmu2*VG/nG+pV2*2*(VS^2)/nS+pV2*2*(VG^2)/nG)
  Cp <- pmuS^2*VS/nS+pmuG^2*VG/nG+pV^2*2*(VS^2)/nS+pV^2*2*(VG^2)/nG
  
  nc <- hp<=p0+Bp+qnorm(1-b)*sqrt(Cp)
  out <- list(hp=hp,Bp=Bp,Cp=Cp,nc=nc,d=d,p0=p0,b=b)
  
  return(out)
}

Senin <- function(muS,muG,VS,VG){
  cvS <- muS/sqrt(VS)
  cvG <- muG/sqrt(VG)
  out <- list(cvS=cvS,cvG=cvG,r=cvS/cvG)
  return(out)
}

Respe <- function(muS,muG,nS,nG,VS,VG,method=1,a=0.05,b=0.2){
  K <- length(nS)
  if (length(nG)!=K){
    muG <- rep(muG[1],K)
    nG <- rep(nG[1],K)
    VG <- rep(VG[1],K)
  }
  
  TT <- (muS-muG)/sqrt(VS/nS+VG/nG)
  if (method==1) {
    ul <- apply( matrix( c( abs(TT-qnorm(1-b/2)),abs(TT+qnorm(1-b/2)) ) ,ncol=2),1,max )
  } else if (method==2) {
    ul <- rep(0,K)
    for (i in 1:K){        
      ul[i] <- uniroot(function(x) pnorm(x-TT[i])-pnorm(-x-TT[i])-1+b, lower=-4-2*abs(TT[i])+TT[i] , upper=4+2*abs(TT[i])+TT[i] )$root
    }
  } else {
    print("Method Type Wrong!")
  }
  ps <- pnorm(qnorm(1-a/2)-ul)-pnorm(-qnorm(1-a/2)-ul)
  out <- list(TT=TT,ul=ul,ps=ps)
  
  return(out)
}

RespeSup <- function(muS,muG,nS,nG,VS,VG,q=0.5,a=0.05,b=0.2){
  K <- length(nS)
  if (length(nG)!=K){
    muG <- rep(muG[1],K)
    nG <- rep(nG[1],K)
    VG <- rep(VG[1],K)
  }
  
  TT <- (muS-q*muG)/sqrt(VS/nS+q^2*VG/nG)
  ul <- TT-qnorm(1-b)
  ps <- pnorm(ul-qnorm(1-a))
  out <- list(TT=TT,ul=ul,ps=ps)
  return(out)
}

Gen <- function(mu1G,mu2G,V1G,V2G,n1G,n2G,mu1,mu2,V1,V2,n1,n2,a=0.05,b=0.2){
  TG <- (mu1G-mu2G)/sqrt(V1G/n1G+V2G/n2G)
  pG <- pnorm(TG-qnorm(1-a/2))-pnorm(-TG-qnorm(1-a/2))
  psG1 <- pnorm(abs(TG)-qnorm(1-b/2)-qnorm(1-a/2))+pnorm(-abs(TG)+qnorm(1-b/2)-qnorm(1-a/2))
  psG <- min(psG1,abs(TG)>qnorm(1-b/2))
  
  TS <- (mu1-mu2)/sqrt(V1/n1+V2/n2)
  pS <- pnorm(TS-qnorm(1-a/2))-pnorm(-TS-qnorm(1-a/2))
  psS1 <- pnorm(abs(TS)-qnorm(1-b/2)-qnorm(1-a/2))+pnorm(-abs(TS)+qnorm(1-b/2)-qnorm(1-a/2))
  psS <- apply( matrix( c(psS1,abs(TS)>qnorm(1-b/2)) ,ncol=2) ,1,min)
  
  out <- list(TG=TG,pG=pG,psG=psG,TS=TS,pS=pS,psS=psS,rp=pS/pG,rps=psS/psG)
  return(out)
}

GenSup <- function(mu1G,mu2G,V1G,V2G,n1G,n2G,mu1,mu2,V1,V2,n1,n2,a=0.05,b=0.2){
  TG <- (mu1G-mu2G)/sqrt(V1G/n1G+V2G/n2G)
  pG <- pnorm(TG-qnorm(1-a))
  psG <- pnorm(TG-qnorm(1-b)-qnorm(1-a))
  
  TS <- (mu1-mu2)/sqrt(V1/n1+V2/n2)
  pS <- pnorm(TS-qnorm(1-a))
  psS <- pnorm(TS-qnorm(1-b)-qnorm(1-a))
  
  out <- list(TG=TG,pG=pG,psG=psG,TS=TS,pS=pS,psS=psS,rp=pS/pG,rps=psS/psG)
  return(out)
}

Gent <- function(mu1G,mu2G,V1G,V2G,n1G,n2G,mu1,mu2,V1,V2,n1,n2,a=0.05,b=0.2){
  K <- length(mu1)
  if (K!=length(mu2)){
    mu2 <- rep(mu2[1],K)
    n2 <- rep(n2[1],K)
    V2 <- rep(V2[1],K)
  }
  
  TG <- (mu1G-mu2G)/sqrt(V1G/n1G+V2G/n2G)
  pG <- 1-pt(qt(1-a/2,n1G+n2G-2),n1G+n2G-2,TG)+pt(-qt(1-a/2,n1G+n2G-2),n1G+n2G-2,TG) 
  llGb <- uniroot(function(x) pt(TG,n1G+n2G-2,x)-1+b/2, upper= 4+TG,lower= -4-4*abs(TG))$root
  ulGb <- uniroot(function(x) pt(TG,n1G+n2G-2,x)-b/2, upper= 4+4*abs(TG),lower= -4+TG)$root
  if ((llGb<=0)*(ulGb>=0)){
    psG <- 0
  } else {
    llG <- min(abs(llGb),abs(ulGb))
    psG <- 1-pt(qt(1-a/2,n1G+n2G-2),n1G+n2G-2,llG)+pt(-qt(1-a/2,n1G+n2G-2),n1G+n2G-2,llG)
  }
  
  TS <- (mu1-mu2)/sqrt(V1/n1+V2/n2)
  pS <- 1-pt(qt(1-a/2,n1+n2-2),n1+n2-2,TS)+pt(-qt(1-a/2,n1+n2-2),n1+n2-2,TS)
  psS <- rep(0,K)
  for (i in 1:K){
    llb <- uniroot(function(x) pt(TS[i],n1[i]+n2[i]-2,x)-1+b/2, upper= 4+TS[i],lower= -4-4*abs(TS[i]))$root
    ulb <- uniroot(function(x) pt(TS[i],n1[i]+n2[i]-2,x)-b/2, upper= 4+4*abs(TS[i]),lower= -4+TS[i])$root
    if ((llb<=0)*(ulb>=0)){
      psS[i] <- 0
    } else {
      ll <- min(abs(llb),abs(ulb))
      psS[i] <- 1-pt(qt(1-a/2,n1+n2-2),n1+n2-2,ll)+pt(-qt(1-a/2,n1+n2-2),n1+n2-2,ll)
    }
  }
  
  out <- list(TG=TG,pG=pG,psG=psG,TS=TS,pS=pS,psS=psS,rp=pS/pG,rps=psS/psG)
  return(out)
}

GentSup <- function(mu1G,mu2G,V1G,V2G,n1G,n2G,mu1,mu2,V1,V2,n1,n2,a=0.05,b=0.2){
  K <- length(mu1)
  if (K!=length(mu2)){
    mu2 <- rep(mu2[1],K)
    n2 <- rep(n2[1],K)
    V2 <- rep(V2[1],K)
  }
  
  TG <- (mu1G-mu2G)/sqrt(V1G/n1G+V2G/n2G)
  pG <- 1-pt(qt(1-a,n1G+n2G-2),n1G+n2G-2,TG)
  llGb <- uniroot(function(x) pt(TG,n1G+n2G-2,x)-1+b, upper= 4+TG,lower= -4-4*abs(TG))$root
  psG <- 1-pt(qt(1-a,n1G+n2G-2),n1G+n2G-2,llGb)
  
  TS <- (mu1-mu2)/sqrt(V1/n1+V2/n2)
  pS <- 1-pt(qt(1-a,n1+n2-2),n1+n2-2,TS)
  llSb <- rep(0,K)
  for (i in 1:K){
    llSb[i] <- uniroot(function(x) pt(TS[i],n1[i]+n2[i]-2,x)-1+b, upper= 4+TS[i],lower= -4-4*abs(TS[i]))$root
  }
  psS <- 1-pt(qt(1-a,n1+n2-2),n1+n2-2,llSb)
  
  out <- list(TG=TG,pG=pG,psG=psG,TS=TS,pS=pS,psS=psS,rp=pS/pG,rps=psS/psG)
  return(out)
}


# ni -> non-inferiority
chow0 <- function(mu0,v0,mud,vd,ni=0){
  r<-seq(0,1,.1)
  pp<-r
  for (i in 1:length(r)){
    md<-r[i]+(1-r[i])*exp(-(mud-mu0)^2/2/(vd+v0))/sqrt(2*pi*(vd+v0))
    a<-(v0+vd)/(v0*vd)
    b<-(-2)*(mu0*vd+mud*v0)/(v0*vd)
    c<-(mu0^2*vd+mud^2*v0)/(v0*vd)
    h<-(-b)/2/a
    k<-c-a*h^2
    p1<-1-pnorm((ni-mud)/sqrt(vd))
    p2<-exp((-1)*k/2)*(1-pnorm((ni-h)*sqrt(a)))/sqrt(a)/sqrt(vd*v0)/sqrt(2*pi)
    pp[i]<-(r[i]*p1+(1-r[i])*p2)/md
  }
  pp # posterior probability
}

# only apply to matched-pair parallel design
chowpair <- function(mu0,v0,muS,VS,nS,ni=0){
  pp <- matrix(0,nrow=K,ncol=11)
  for (i in 1:K){
    pp[i,] <- chow0(mu0[i],v0[i],muS[i],VS[i]/nS[i],ni)
  }
  pp
}


Japan <- function(muS,muG,rho=0.8){
  nc <- 1 - ( (muS*muG)>=0 )*( abs(muS)>=(abs(muG)*rho) )
  return(nc)
}


CIllmp <- function(mu,V,a,n,df,side=1){
  if (side==1){
    aa <- a
  }else if (side==2){
    aa <- a/side
  }else{
    print("Wrong side!")
  }
  cill <- mu-sqrt(V/n)*qt(1-aa,df)
  return <- cill
}

CIll2gVeq <- function(muT,muR,VT,VR,nT,nR,a,side=1){
  if (side==1){
    aa <- a
  }else if (side==2){
    aa <- a/side
  }else{
    print("Wrong side!")
  }
  cill <- (muT-muR)-sqrt( (VT*(nT-1)+VR*(nR-1))/(nT+nR-2)*((nT+nR)/nT/nR) )*qt(1-aa,nT+nR-2)
  return <- cill
}

# categmean <- function(mat,categ,K,nSdif,m){
# 	ncateg <- matrix(0,nrow=1+m,ncol=1+nSdif)
# 	ncateg[1,] <- c(nrow(mat)*(K-nSdif),rep(nrow(mat),nSdif))
# 	mat0 <- mat[,1:(K-nSdif)]
# 	categ0 <- categ[,1:(K-nSdif)]
# 	mat1 <- mat[,(K-nSdif+1):K]
# 	categ1 <- categ[,(K-nSdif+1):K]
# 	vec0 <- as.vector(mat0)
# 	mean0 <- c(mean(vec0),apply(mat[,(K-nSdif+1):K],2,mean))
# 	mean0 <- matrix(mean0,nrow=1)
# 	mean1 <- matrix(0,nrow=m,ncol=nSdif+1)
# 	ncategall <- matrix(0,nrow=1+m,ncol=1)
# 	ncategall[1] <- nrow(mat)*K
# 	meanall <- matrix(0,nrow=1+m,ncol=1)
# 	meanall[1] <- mean(mat)
# 	for (i in 1:m){
# 		ncategall[i+1] <- sum(categ==i)
# 		meanall[i+1] <- mean(mat[categ==i])
# 		ncateg[i+1,1] <- sum(categ0==i)
# 		ncateg[i+1,2:(1+nSdif)] <- apply(categ1==i,2,sum)
# 		mean1[i,1] <- mean(mat0[categ0==i])
# 		for (j in 1:nSdif){
# 			mean1[i,1+j] <- mean(mat1[categ1[,j]==i,j])
# 		}
# 	}
# 	mean2 <- rbind(mean0,mean1)
# 	mean3 <- cbind(matrix(seq(0,m),ncol=1),round(mean2,3),round(meanall,3),ncateg,ncategall)
# 	return(mean3)
# }

categmean <- function(mat,categ,cK,m){
  ind2 <- cumsum(cK)
  ind1 <- ind2+1-cK
  C <- length(cK)
  ncateg <- matrix(0,nrow=1+m,ncol=C+1)
  ncategp <- matrix(1,nrow=1+m,ncol=C+1)
  meanall <- matrix(0,nrow=1+m,ncol=C+1)
  for (i in 1:m){
    ncateg[i+1,C+1] <- sum(as.vector(categ)==i)
    meanall[i+1,C+1] <- mean(mat[categ==i])
  }
  for (j in 1:C){
    ncateg[1,j] <- sum(as.vector(categ[,ind1[j]:ind2[j]])!=0)
    tmp1 <- as.vector(mat[,ind1[j]:ind2[j]])
    tmp2 <- tmp1[as.vector(categ[,ind1[j]:ind2[j]])!=0]
    meanall[1,j] <- mean(tmp2)
    for (i in 1:m){
      ncateg[1+i,j] <- sum(as.vector(categ[,ind1[j]:ind2[j]])==i)
      tmp1 <- as.vector(mat[,ind1[j]:ind2[j]])
      tmp2 <- tmp1[as.vector(categ[,ind1[j]:ind2[j]])==i]
      meanall[1+i,j] <- mean(tmp2)
    }
  }
  ncateg[1,C+1] <- sum(ncateg[1,1:C])
  meanall[1,C+1] <- mean(mat[categ!=0])
  ncategp[2:(1+m),] <- ncateg[2:(1+m),]/matrix(ncateg[1,],nrow=m,ncol=C+1,byrow=TRUE)
  
  head <- matrix(seq(0,m),nrow=1+m,ncol=1)
  mean3 <- cbind(head,meanall,ncateg,ncategp)
  return(mean3)
}
