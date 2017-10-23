### Survival Analysis - exponential Model ###

## Input Parameters ##

# sample size calculation (N-query/Clinical protocol)

t.actual <- 2
t.followup <- 1
t.mean <- 0.5*t.actual+t.followup

refer.event.rate <- 0.02
alpha <- 0.05
beta <- 0.1
power <- 1-beta

relative.risk.reduc <- 0.2

#the sample size under power <-  1-beta <-  0.9
N <-19520
n <- 0.2*N

# exponential model parameters (constant hazard model):
lambda.g.trt <- 0.1
lambda.g.ctr <- 0.15

lambda.g.trt = 0.0081
lambda.g.ctr = 0.0101
# Here we already assume the constant hazard in treatment group is higher than the control group.

## Simulation ##

#generate random lognormal sample
seed=111
set.seed(seed)
global.trt <- rexp(N/2,rate = lambda.g.trt)
global.ctr <- rexp(N/2,rate = lambda.g.ctr)
lifetime.g <- c(global.trt,global.ctr)
grp.g <- c(rep(1,N/2),rep(0,N/2))
   # 1 for treatment, 0 for control
delta.g <- as.numeric(lifetime.g<=t.mean)
   # 1 for failure, 0 for censored
global <- data.frame(lifetime.g,delta.g,grp.g) 
   # combine survival data together.

set.seed(seed)
trt <- global[which(global$grp.g==1),]
ctr <- global[which(global$grp.g==0),]
# keep 1 to 1 ratio between trt and ctr in region
region.trt <- trt[sample(nrow(trt),n/2,replace = FALSE),]
region.ctr <- ctr[sample(nrow(ctr),n/2,replace = FALSE),]
region <- rbind(region.trt,region.ctr)



library(survival)
# log-rank test
  # test global significance
logrank <- survdiff(Surv(lifetime.g,delta.g) ~ grp.g)
summary(logrank)
# p= 0 No signidicance difference in survival curve.

# cox-proportional hazard model
cph.g <- coxph(Surv(lifetime.g,delta.g)~ grp.g,data = global)
cph.r <- coxph(Surv(lifetime.g,delta.g)~ grp.g,data = region)

sum.cph.g <- summary(cph.g)
sum.cph.r <- summary(cph.r)


hard.g <- exp(cph.g$coefficients)
hard.r <- exp(cph.r$coefficients)

sd.g <- sqrt(cph.g$var[1,1])
sd.r <- sqrt(cph.r$var[1,1])

upper.g <- exp(cph.g$coefficients+qnorm(1-alpha/2)*sd.g)
upper.r <- exp(cph.r$coefficients+qnorm(1-alpha/2)*sd.r)

#superiority test
#presumption: hard.g<1
if(hard.r<hard.g){
   flag=1
}else if((hard.g<hard.r<1)&(upper.r<1)){
   flag=2
}else if((hard.g<hard.r<1)&(upper.r>1)){
   flag=3
}else if(hard.g<1<=hard.r){
   flag=4
}


