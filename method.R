## same consistency method
## different lambda


## read in prop.css data

# lambda=1

# lambda=2

# lambda=3
doc <- "prop.css.3.2017-07-23"
file.address <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/lambda/",
                      doc,
                      ".txt",
                      sep = "")

lambda.3 <- read.table(file=file.address,
                 header = FALSE,
                 sep = " ",
                 skip=1)
(lambda.3  <- lambda.3[,2:4])
(colnames(lambda.3) <- c('Bayesian.min','Bayesian.mean','Japanese'))
lambda.3

# lambda=4
# lambda=5

doc <- "prop.css.5.2017-07-22"
file.address <- paste("C:/Users/jings/Dropbox/Internship/BMS/Intern project/lambda/",
                      doc,
                      ".txt",
                      sep = "")

lambda.5 <- read.table(file=file.address,
                       header = FALSE,
                       sep = " ",
                       skip=1)
(lambda.5  <- lambda.5[,2:4])
(colnames(lambda.5) <- c('Bayesian.min','Bayesian.mean','Japanese'))
lambda.5


## Plot!

library(gcookbook)
library(ggplot2)
n.percent <- length(percent.of.g)

## Bayesian.min
m<-1

lambda.grp <- c(rep('0.01',n.percent),
                 rep('0.05',n.percent),
                 rep('0.1',n.percent),
                 rep('0.5',n.percent),
                 rep('0.8',n.percent))

percent.css <- rep(percent.of.g,5)
prop.css.vector <- c(lambda.1[,m],
                     lambda.2[,m],
                     lambda.3[,m],
                     lambda.4[,m],
                     lambda.5[,m])
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

css +  guides(color=guide_legend(title="Method")) 


## Bayesian.mean
m<-2

## Japanese
m<-3
