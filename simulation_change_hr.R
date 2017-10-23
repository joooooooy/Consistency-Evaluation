# fixed lambda.ctr = 0.1
# fixed pwr=0.9
# change HR:
# we already test for HR=0.8 
k=1
k=3
hr <- c(0.7,0.8,0.9)
N <- 2*c(1084,2590,10980)
time.hr <- rep(NA,length(hr))


## Nquery Input/Output
## Log-rank test of survival in two groups followed for fixed time

lambda.g.trt <- c(0.07,0.08,0.09)
lambda.g.ctr <- c(0.01,0.05,0.1,0.5,0.8) #we only use lambda.ctr=0.1 here

N[k]
(percent.of.g <- seq(from=0.01,to=0.99,by=0.01))
(sample.size <- percent.of.g*N[k])


seed.from <- 1 
seed.to <-20000
seed <- seed.from:seed.to 

## simulation under pwr[k]

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
# value of jp2 : 0 1
# consistent 1 for rrr.r>=rrr.g*0.5


##################################
######## TIME ESTIMATION #########
##################################

i=1
t1<-Sys.time()
for(J in 1:10){
  temp <- consistency(N=N[k],
                      n=sample.size[i],
                      seed=seed[j],
                      lambda.g.trt = lambda.g.trt[k],
                      lambda.g.ctr = lambda.g.ctr[3])
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
t2<-Sys.time()
t2-t1
20000*99*(t2-t1)/10
20000*99*(t2-t1)/10/60/60


##################################
##############  END ##############
##################################
t1<-Sys.time()
for(i in 1:length(percent.of.g)){
  for(j in 1:length(seed)){
    temp <- consistency(N=N[k],
                        n=sample.size[i],
                        seed=seed[j],
                        lambda.g.trt = lambda.g.trt[k],
                        lambda.g.ctr = lambda.g.ctr[3])
    HR.class[i,j] <- temp$flag
    if(temp$upper.g<1){
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
(time.hr[k]=t2-t1)


######################################
### ??????simulation???????????????????????? ###
######################################


## Store simulation results on cvs files
(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "HR.class.hr=",as.character(hr[k]),".",Sys.Date(),".csv",sep = ""))
write.csv(x=HR.class,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               " byes.min.hr=",as.character(hr[k]),".",Sys.Date(),".csv",sep = ""))
write.csv(x= byes.min,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "byes.mean.hr=",as.character(hr[k]),".",Sys.Date(),".csv",sep = ""))
write.csv(x=byes.mean,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "jp.hr=",as.character(hr[k]),".",Sys.Date(),".csv",sep = ""))
write.csv(x=jp,file = name)

## HR categories 1 2 3 4
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
thred <- 0.8
for(i in 1:length(percent.of.g)){
  for(j in 1:length(seed)){
    byes.min.cat[i,j] <- as.numeric(byes.min[i,j] >= thred)
  }
}
freq.count.byes.min <- t(apply(X = byes.min.cat,MARGIN = 1,FUN = table))
freq.count.byes.min

for(i in 1:length(percent.of.g)){
  for(j in 1:length(seed)){
    byes.mean.cat[i,j] <- as.numeric(byes.mean[i,j] >= thred)
  }
}
freq.count.byes.mean <- t(apply(X = byes.mean.cat,MARGIN = 1,FUN = table))
freq.count.byes.mean

## Japanese

freq.count.jp <- t(apply(X = jp,MARGIN = 1,FUN = table))
freq.count.jp


## store frequency information
(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "freq.count.hr.hr=",as.character(hr[k]),".",Sys.Date(),".txt",sep = ""))
write.table(freq.count.hr,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "freq.count.byes.min.hr=",as.character(hr[k]),".",Sys.Date(),".txt",sep = ""))
write.table(freq.count.byes.min,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "freq.count.byes.mean.hr=",as.character(hr[k]),".",Sys.Date(),".txt",sep = ""))
write.table(freq.count.byes.mean,file = name)

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "freq.count.jp.hr=",as.character(hr[k]),".",Sys.Date(),".txt",sep = ""))
write.table(freq.count.jp,file = name)

# the proportion of each classification of the result under varying power.

for(i in 1:length(percent.of.g)){
  prop.hr[i,]<- freq.count.hr[i,]/(length(seed)-freq.count.hr[i,4])
}

(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "prop.hr.hr=",as.character(hr[k]),".",Sys.Date(),".txt",sep = ""))
write.table(x=prop.hr,file = name)


## ploting.hr
#install.packages('gcookbook')
#install.packages('ggplot2')
library(gcookbook)
library(ggplot2)

n.percent <- length(percent.of.g)
hard.grp <- c(rep('1',n.percent),
              rep('3',n.percent),
              rep('4',n.percent),
              #rep('999',n.percent),
              rep('2',n.percent))
(percent.of.g <- seq(from=0.01,to=0.99,by=0.01))
percent.hr <- rep(percent.of.g,4)
prop.hr.vector <- c(prop.hr[,1],
                    prop.hr[,2],
                    prop.hr[,3],
                    #prop.hr[,4],
                    prop.hr[,5])
(combined.hr <- data.frame(hard.grp,
                           percent.hr,
                           prop.hr.vector))

note <- paste("Control Group ~ EXP(0.1), power=0.8, HR=",hr[k],sep="")

hr <- ggplot(combined.hr, 
             aes(x=factor(percent.hr), 
                 y=prop.hr.vector, 
                 colour=hard.grp,
                 group=hard.grp)) + geom_line(size=2) 
hr <- hr + labs(x = "Regional Sample Size",y="Percentage",
                title="Hazard Ratio Categories",
                subtitle="Consistency Evaluation",
                caption=note)
hr <- hr+ scale_x_discrete(breaks=seq(0,1,0.05)) 

hr +  guides(color=guide_legend(title=NULL)) 


## merge the result of bayes and japan1&2

hr <- c(0.7,0.8,0.9)
(name <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/",
               "prop.css.hr=",as.character(hr[k]),".",Sys.Date(),".txt",sep = ""))
write.table(x=prop.css,file = name)

## ploting.byes.jp
n.percent <- length(percent.of.g)
consist.grp <- c(rep('Bayesian.min',n.percent),
                 rep('Bayesian.mean',n.percent),
                 rep('japanese',n.percent))

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

css <- css + labs(x = "Regional Sample Size",y="Percentage",
                  title="Bayesian Approach Compared to Japanese Approach",
                  subtitle="Consistency Evaluation",
                  caption=note)

css <- css+ scale_x_discrete(breaks=seq(0,1,0.05))  

css <- css + guides(color=guide_legend(title="Method")) 

css <- css + expand_limits(y=0)

css <- css + scale_y_continuous(breaks = seq(0,1,0.2))

css

install.packages("scales")
library(scales)

