## Nquery Input/Output
## Log-rank test of survival in two groups followed for fixed time

## Scenerio 1 
## fixed HR = 0.8 power=90%
## change exponential parameter in control group 0.01,0.05,0.1,0.5,0.8
lambda.g.trt <- c(0.008,0.04,0.08,0.4,0.64)
lambda.g.ctr <- c(0.01,0.05,0.1,0.5,0.8)
(N <- 2*c(23621,4956,2590,720,560))
num.event.req <- rep(844,5)
time.lambda <- rep(NA,length(lambda.g.ctr))


## k=1,2,3,4,5
k=1
N[k]
(percent.of.g <- seq(from=0.01,to=0.99,by=0.01))
(sample.size <- percent.of.g*N[k])

# set seed
seed.from <- 1 
seed.to <-20000
seed <- seed.from:seed.to

## simulation under lambda.g.ctr[k]

# The result of each simulation
HR.class <- matrix(data = NA,nrow = length(percent.of.g),ncol= length(seed))
# value of HR.class : 1 2 3 4 999
byes.min <- matrix(data = NA,nrow = length(percent.of.g),ncol= length(seed))
byes.min.cat <- matrix(data = NA,nrow = length(percent.of.g),ncol= length(seed))
# value of byes : 0 1
# consistent 1 for (min(pp)>=0.8)  
byes.mean <- matrix(data = NA,nrow = length(percent.of.g),ncol= length(seed))
byes.mean.cat <- matrix(data = NA,nrow = length(percent.of.g),ncol= length(seed))
# value of byes : 0 1
# consistent 1 for (mean(pp)>=0.8) 
jp <- matrix(data = NA,nrow = length(percent.of.g),ncol= length(seed))
# value  of jp2 : 0 1
# consistent 1 for rrr.r>=rrr.g*rho 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # testing consistency() before simulation loop # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
(temp <- consistency(N=N[k],
                     n=sample.size[55],
                     seed=111,
                     lambda.g.trt = lambda.g.trt[k],
                     lambda.g.ctr = lambda.g.ctr[k]))
# notice here that flag =1 (hard.r<hard.g)
# upper.g < 1 while upper.r>1
min(chow0(mu0 = temp$u.g,
      v0 = temp$var.g,
      mud = temp$u.r,
      vd = temp$var.r,
      ni=0))

mean(chow0(mu0 = temp$u.g,
          v0 = temp$var.g,
          mud = temp$u.r,
          vd = temp$var.r,
          ni=0))

Japan2(hard.g = temp$hard.g,
       hard.r = temp$hard.r,
       rho = 0.5)


i=1
t1<-Sys.time()
for(J in 1:10){
  temp <- consistency(N=N[k],
                      n=sample.size[i],
                      seed=seed[j],
                      lambda.g.trt = lambda.g.trt[k],
                      lambda.g.ctr = lambda.g.ctr[k])
  HR.class[i,j] <- temp$flag
  
  ##############
  ####UPDATE####
  ##############
  if(temp$upper.g<1){
    
    # when global significant
    byes.pp <- chow0(mu0 = temp$u.g,
                     v0 = temp$var.g,
                     mud = temp$u.r,
                     vd = temp$var.r,
                     ni=0)
    byes.min[i,j] <- min(byes.pp)
    byes.mean[i,j] <- mean(byes.pp)
    jp[i,j] <- Japan2(hard.g = temp$hard.g,
                      hard.r = temp$hard.r,
                      rho = 0.5)
  }
}
t2<-Sys.time()
t2-t1
20000*99*(t2-t1)/10
20000*99*(t2-t1)/10/60/60
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# expect to take 100 hr
# end: "2017-07-29 11:23:00 EDT"

t1<-Sys.time()
for(i in 55:length(percent.of.g)){
  for(j in 1:length(seed)){
    temp <- consistency(N=N[k],
                        n=sample.size[i],
                        seed=seed[j],
                        lambda.g.trt = lambda.g.trt[k],
                        lambda.g.ctr = lambda.g.ctr[k])
    HR.class[i,j] <- temp$flag

    if(temp$upper.g<1){
      
      # when global significant
      byes.pp <- chow0(mu0 = temp$u.g,
                       v0 = temp$var.g,
                       mud = temp$u.r,
                       vd = temp$var.r,
                       ni=0)
      byes.min[i,j] <- min(byes.pp)
      byes.mean[i,j] <- mean(byes.pp)
      jp[i,j] <- Japan2(hard.g = temp$hard.g,
                        hard.r = temp$hard.r,
                        rho = 0.5)
    }


  }
}
t2<-Sys.time()
(time.lambda[k]=t2-t1)

############################################################################
##  ###  ###  ###  ###  ###  ###  ### ###  ###  ###  ###  ###  ###  ###  ###
############################################################################


## Store simulation results on cvs files
(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "HR.class.",as.character(k),".",Sys.Date(),".csv",sep = ""))
write.csv(x=HR.class,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "byes.min.",as.character(k),".",Sys.Date(),".csv",sep = ""))
write.csv(x= byes.min,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "byes.mean.",as.character(k),".",Sys.Date(),".csv",sep = ""))
write.csv(x=byes.mean,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "jp.",as.character(k),".",Sys.Date(),".csv",sep = ""))
write.csv(x=jp,file = name)

## HR categories 1 2 3 4 999
## calculate frequency
freq.count.hr <- apply(X = HR.class,MARGIN = 1,FUN = table)
freq.count.hr

library(gtools)
temp.freq.hr <- freq.count.hr[[1]]
for(i in 2:length(percent.of.g)){
  temp.freq.hr <- smartbind(temp.freq.hr,freq.count.hr[[i]])
}
freq.count.hr<-temp.freq.hr
row.names(freq.count.hr) <- as.character(1:length(percent.of.g))
freq.count.hr[is.na(freq.count.hr)]<-0
freq.count.hr



## bayesian
## calculate frequency
thred <- 0.8

for(i in 1:length(percent.of.g)){
  for(j in 1:length(seed)){
    byes.min.cat[i,j] <- as.numeric(byes.min[i,j] >= thred)
  }
}
freq.count.byes.min <- t(apply(X = byes.min.cat,MARGIN = 1,FUN = table))
freq.count.byes.min
row.names(freq.count.byes.min) <- as.character(1:length(percent.of.g))
freq.count.byes.min
##############
####UPDATE####
##############
freq.count.byes.min <- apply(X = byes.min.cat,MARGIN = 1,FUN = table)
temp.freq.byes.min <- freq.count.byes.min[[1]]
for(i in 2:length(percent.of.g)){
  temp.freq.byes.min <- smartbind(temp.freq.byes.min,freq.count.byes.min[[i]])
}
(freq.count.byes.min<-temp.freq.byes.min)
row.names(freq.count.byes.min) <- as.character(1:length(percent.of.g))
freq.count.byes.min[is.na(freq.count.byes.min)]<-0
freq.count.byes.min



for(i in 1:length(percent.of.g)){
  for(j in 1:length(seed)){
    byes.mean.cat[i,j] <- as.numeric(byes.mean[i,j] >= thred)
  }
}
freq.count.byes.mean <- apply(X = byes.mean.cat,MARGIN = 1,FUN = table)
freq.count.byes.mean

##############
####UPDATE####
##############
(temp.freq.byes.mean <- freq.count.byes.mean[[1]])
for(i in 2:length(percent.of.g)){
  temp.freq.byes.mean <- smartbind(temp.freq.byes.mean,freq.count.byes.mean[[i]])
}
(freq.count.byes.mean<-temp.freq.byes.mean)
row.names(freq.count.byes.mean) <- as.character(1:length(percent.of.g))
freq.count.byes.mean[is.na(freq.count.byes.mean)]<-0

freq.count.byes.mean



## Japanese
## calculate frequency
freq.count.jp <- t(apply(X = jp,MARGIN = 1,FUN = table))
freq.count.jp
row.names(freq.count.jp) <- as.character(1:length(percent.of.g))
freq.count.jp
##############
####UPDATE####
##############
freq.count.jp <- apply(X = jp,MARGIN = 1,FUN = table)
temp.freq.jp <- freq.count.jp[[1]]
for(i in 2:length(percent.of.g)){
  temp.freq.jp <- smartbind(temp.freq.jp,freq.count.jp[[i]])
}
(freq.count.jp<-temp.freq.jp)
row.names(freq.count.jp) <- as.character(1:length(percent.of.g))
freq.count.jp[is.na(freq.count.jp)]<-0
freq.count.jp



## store frequency information
(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "freq.count.hr.",as.character(k),".",Sys.Date(),".txt",sep = ""))
write.table(freq.count.hr,file = name)


(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "freq.count.byes.min.",as.character(k),".",Sys.Date(),".txt",sep = ""))
write.table(freq.count.byes.min,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "freq.count.byes.mean.",as.character(k),".",Sys.Date(),".txt",sep = ""))
write.table(freq.count.byes.mean,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "freq.count.jp.",as.character(k),".",Sys.Date(),".txt",sep = ""))
write.table(freq.count.jp,file = name)

# the proportion of each classification of the result under varying power.
prop.hr <- matrix(data = NA, nrow = length(percent.of.g),ncol = 5)
for(i in 1:length(percent.of.g)){
  prop.hr[i,]<- as.numeric(freq.count.hr[i,]/(length(seed)-freq.count.hr[i,4]))
}

row.names(prop.hr) <- as.character(1:99)
colnames(prop.hr) <- c("1","3","4","999","2")
prop.hr

prop.hr[91:97,3]<-0

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "prop.hr.",as.character(k),".",Sys.Date(),".txt",sep = ""))
write.table(x=prop.hr,file = name)


## ploting.hr
#install.packages('gcookbook')
#install.packages('ggplot2')
library(gcookbook)
library(ggplot2)

n.percent <- length(percent.of.g)-3
hard.grp <- c(rep('1',n.percent),
              rep('3',n.percent),
              rep('4',n.percent),
              #rep('999',n.percent),
              rep('2',n.percent))
(percent.of.g <- seq(from=0.01,to=0.96,by=0.01))
percent.hr <- rep(percent.of.g,4)
prop.hr.vector <- c(prop.hr[1:n.percent,1],
                    prop.hr[1:n.percent,2],
                    prop.hr[1:n.percent,3],
                    #prop.hr[,4],
                    prop.hr[1:n.percent,5])
(combined.hr <- data.frame(hard.grp,
                           percent.hr,
                           prop.hr.vector))

combined.hr[289:298,3]<-0

note <- paste("Control Group ~ EXP(",lambda.g.ctr[k],")"," power=0.9, HR=0.8",sep = "")
note

hr <- ggplot(combined.hr, 
       aes(x=factor(percent.hr), 
           y=prop.hr.vector, 
           colour=hard.grp,
           group=hard.grp)) + geom_line(size=2) 
hr <- hr + labs(x = "Proportion of Regional Sample Size",y="Proportion of HR category",
        title="Hazard Ratio Categories",
        subtitle="Consistency Evaluation",
        caption=note)
hr <- hr+ scale_x_discrete(breaks=seq(0,1,0.05)) 

hr +  guides(color=guide_legend(title=NULL)) 


##############
####UPDATE####
##############

## merge the result of bayes and japan1&2

## prop.byes.min <- freq.count.byes.min[,2]/length(seed)
prop.byes.min <- matrix(data = NA,nrow = length(percent.of.g),ncol = 2)
for(i in 1:length(percent.of.g)){
  prop.byes.min[i,]<- as.numeric(freq.count.byes.min[i,]/sum(freq.count.byes.min[i,]))
}
row.names(prop.byes.min) <- as.character(1:length(percent.of.g))
colnames(prop.byes.min) <- as.character(0:1)
prop.byes.min

## prop.byes.mean <- freq.count.byes.mean[,2]/length(seed)
prop.byes.mean <- matrix(data = NA,nrow = length(percent.of.g),ncol = 2)
for(i in 1:length(percent.of.g)){
  prop.byes.mean[i,]<-as.numeric(freq.count.byes.mean[i,]/sum(freq.count.byes.mean[i,]))
}
row.names(prop.byes.mean) <- as.character(1:length(percent.of.g))
colnames(prop.byes.mean) <- as.character(0:1)
prop.byes.mean

#prop.jp <- freq.count.jp[,2]/length(seed)
prop.jp <- matrix(data = NA,nrow = length(percent.of.g),ncol = 2)
for(i in 1:length(percent.of.g)){
  prop.jp[i,]<-as.numeric(freq.count.jp[i,]/sum(freq.count.jp[i,]))
}
row.names(prop.jp) <- as.character(1:length(percent.of.g))
colnames(prop.jp) <- as.character(0:1)
prop.jp


(prop.css <- cbind(prop.byes.min[,2],prop.byes.mean[,2],prop.jp[,2]))

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "prop.css.",as.character(k),".",Sys.Date(),".txt",sep = ""))
write.table(x=prop.css,file = name)

## ploting.byes.jp1.jp2
n.percent <- length(percent.of.g)
consist.grp <- c(rep('Bayesian.min',n.percent),
                 rep('Bayesian.mean',n.percent),
                 rep('Japanese',n.percent))

percent.css <- rep(percent.of.g,3)
prop.css.vector <- c(prop.css[,1],
                     prop.css[,2],
                     prop.css[,3])
(combined.css <- data.frame(consist.grp,
                            percent.css,
                            prop.css.vector))

css <- ggplot(combined.css, 
       aes(x=factor(percent.css), 
           y=prop.css.vector, 
           colour=consist.grp,
           group=consist.grp)) + geom_line(size=2)

css <- css + labs(x = "Proportion of Regional Sample Size",y="Proportion of Consistent Result",
                title="Bayesian Approach Compared to Japanese Approach",
                subtitle="Consistency Evaluation",
                caption=note)
css <- css+ scale_x_discrete(breaks=seq(0,1,0.05)) 

css +  guides(color=guide_legend(title="Method")) 
