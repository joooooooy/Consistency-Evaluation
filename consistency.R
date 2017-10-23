# N - the global sample size
# n - the regional sample size
# alpha - Significance level
# lambda.g.trt - exponential parameter for global treatment group
# lambda.g.ctr - exponential parameter for global control group
consistency <- function(N,n,seed,
                        alpha=0.05,
                        lambda.g.trt,lambda.g.ctr){
  t.actual <- 2
  t.followup <- 1
  t.mean <- 0.5*t.actual+t.followup
  
  #set.seed(seed)
  #??# previous simulation result is based on seed <- 1:1000
  #??# I am not sure if it somehow leads to the same pattern
  #??# Now I am testing the simulation without seed now
  
  global.trt <- rexp(N/2,rate = lambda.g.trt)
  global.ctr <- rexp(N/2,rate = lambda.g.ctr)
  lifetime.g <- c(global.trt,global.ctr)
  grp.g <- c(rep(1,N/2),rep(0,N/2))
  # 1 for treatment, 0 for control
  delta.g <- as.numeric(lifetime.g<=t.mean)
  # 1 for failure, 0 for censored
  global <- data.frame(lifetime.g,delta.g,grp.g) 
  # combine survival data together.
  num.failure <- sum(delta.g)
  # total number of failure
  
  
  ##############
  ####UPDATE####
  ##############
  trt <- global[1:N/2,]
  ctr <- global[(N/2+1):N,]
  # keep 1 to 1 ratio between trt and ctr in region
  region.trt <- trt[sample(nrow(trt),round(n/2),replace = FALSE),]
  region.ctr <- ctr[sample(nrow(ctr),round(n/2),replace = FALSE),]
  region <- rbind(region.trt,region.ctr)
  
  library(survival)
  # two-sided log-rank test
  # test global significance
  # global.logrank <- survdiff(Surv(lifetime.g,delta.g) ~ grp.g)
  
  # > survdiff(Surv(lifetime.g,delta.g) ~ grp.g)
  #   Call:
  #   survdiff(formula = Surv(lifetime.g, delta.g) ~ grp.g)
  
  # N Observed Expected (O-E)^2/E (O-E)^2/V
  # grp.g=0 23220      465      410      7.26      14.5
  # grp.g=1 23220      357      412      7.24      14.5
  
  # Chisq= 14.5  on 1 degrees of freedom, p= 0.000141 
  
  
  # p.val <- pchisq(q=global.logrank$chisq,df=1,lower.tail = FALSE)
  # p= 0 No signidicance difference in survival curve.
  
  
  # cox-proportional hazard model
  cph.g <- coxph(Surv(lifetime.g,delta.g)~ grp.g,data = global)
  cph.r <- coxph(Surv(lifetime.g,delta.g)~ grp.g,data = region)
  summary(cph.g)
  summary(cph.r)
  
  hard.g <- exp(cph.g$coefficients)
  # exp(coef)
  hard.r <- exp(cph.r$coefficients)
  
  sd.g <- sqrt(cph.g$var[1,1])
  # se(coef)
  sd.r <- sqrt(cph.r$var[1,1])
  
  upper.g <- exp(cph.g$coefficients+qnorm(1-alpha/2)*sd.g)
  upper.r <- exp(cph.r$coefficients+qnorm(1-alpha/2)*sd.r)
  
  
  # verifying lisa's code
  # > cph.g$coefficients
  # grp.g 
  # -0.2670747 
  # > ((-log(upper.g)+cph.g$coefficients)/qnorm(1-alpha/2))^2
  #  grp.g 
  # 0.00495166 
  #  > cph.g$var[1,1]
  # [1] 0.00495166
  
  
  #superiority test
  #presumption: hard.g<1
  if((hard.r<hard.g)&(upper.g<1)){
    flag=1
  }else if((hard.g<hard.r)&(hard.r<1)&(upper.r<1)&(upper.g<1)){
    flag=2
  }else if((hard.g<hard.r)&(hard.r<1)&(upper.r>1)&(upper.g<1)){
    flag=3
  }else if((hard.g<1)&(1<=hard.r)&(upper.g<1)){
    flag=4
  }else{
    flag=999
    }
  
  return(list(flag=flag,
              hard.g=hard.g,
              hard.r=hard.r,
              u.g = -cph.g$coefficients,
              u.r = -cph.r$coefficients,
              var.g = cph.g$var[1,1],
              var.r = cph.r$var[1,1],
              upper.g = upper.g,
              upper.r = upper.r,
              num.event=sum(delta.g)))
}


# ni -> non-inferiority margin
# _0 for global
# _d for regional
# mu for mean 
# v. for variance
# Bayesian approach
chow0 <- function(mu0,v0,mud,vd,ni=0){
  r<-seq(0,1,.1)
  # the weight r reflects relative confidence of the regulatory authority on the evidence provided by
  # the original region would be borrowed.
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
  pp 
  # posterior probability
}


# For normal situation
Japan1 <- function(muG,muS,rho=0.8){
  nc <- 1 - ( (muS*muG)>=0 )*( abs(muS)>=(abs(muG)*rho) )
  return(nc)
}

# For survival situation
Japan2 <- function(hard.g,hard.r,rho=0.5){
  rrr.g <- 1-hard.g
  rrr.r <- 1-hard.r
  consist <- as.numeric(rrr.r>=rrr.g*rho)
  return(consist)
}