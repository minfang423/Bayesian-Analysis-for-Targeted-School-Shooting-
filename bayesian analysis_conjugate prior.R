#-------------------
# load data
#-------------------
library(readxl)
library(moments)
setwd("~/Documents/Min/Study/Math/Statistics/SJSU MATH 264 Bayesian Analysis/project/264 project/school shooting")
data <- read_excel("data.xlsx",sheet = 2)
unique(data$shooting_type)

# create dummy variable for shooting type
data$type = ifelse(data$shooting_type == "targeted" |data$shooting_type == "targeted and indiscriminate",1,0 )
unique(data$type)

#-------------------
# summary
#-------------------
sum(data$type)
# 136
mean(data$type)
# 0.6181818
table(data$type)
#         0         1 
# 0.3818182 0.6181818 

hist(data$type,col = "cornflowerblue",border = "wheat",main="school shooting by type")
text(0.05, 90,"non-targetted",cex = 0.8, col = "black")
text(0.93, 140,"targetted",cex = 0.8, col = "black")
title(xlab = "")


#-------------------
# predictive
# 
# Conjugate Prior: Beta(37.1,16.96)
# sampling: 
# posterior:beta(37.1+sum(y), n-sum(y)+16.96) = beta(173.1,100.96); mean=0.68,mode=0.632581
# predict:bin(1,1,theta)

# #-------------------
library("LearnBayes")
source('~/Documents/Min/Study/Math/Statistics/SJSU MATH 264 Bayesian Analysis/project/264 project/findBeta.R')
quantile1 <- list(p=0.5, x=0.6)    # we believe the median of the prior is 0.85
quantile2 <- list(p=0.99999,x=0.9) # we believe the 99.999th percentile of the prior is 0.95
quantile3 <- list(p=0.00001,x=0.4) # we believe the 0.001st percentile of the prior is 0.60
findBeta(quantile1,quantile2,quantile3)
# [1] "The best beta prior has a= 37.0996896896897 b= 16.9612312312312"



##########################################  1.Frequest CI ########################################## 
library(binom)
binom.confint(sum(data$type), length(data$type), conf.level = 0.95, method = 'all')
#    x   n      mean     lower     upper
#  136 220 0.6181818 0.5524446 0.6798627
# plot
fre_l = 0.5524446
fre_r = 0.6798627
xf = seq(0.4,0.8,length=1000)
y_bar = mean(data$type)
s.sq = y_bar * (1-y_bar)/220
yf = dnorm(xf,mean=y_bar,sd=sqrt(s.sq))
plot(xf,yf,type="l",xlab = "theta",ylab="",main="95% Confidence Interval for theta")

# Returns a vector of boolean values representing whether the x value is between the two bounds then
# filters the values so that only the ones within the bounds are returned
bounds_filter <- xf >= fre_l & xf <= fre_r
xf_within_bounds <- xf[bounds_filter]
yf_within_bounds <- yf[bounds_filter]
# We want the filled in area to extend all the way down to the y axis which is why these two lines are necessary
# It makes the first point in the polygon (lower_bound, 0) and the last point (upper_bound, 0)
xf_polygon <- c(fre_l, xf_within_bounds, fre_r)
yf_polygon <- c(0, yf_within_bounds, 0)
polygon(xf_polygon, yf_polygon, col = "skyblue")
text(fre_l, 4, "0.552", srt = 90, cex = 0.8, col = "gray30")
text(fre_r, 4, "0.679", srt = 90, cex = 0.8, col = "gray30")


########################################## 2. HPD ########################################## 
set.seed(123)
n.sim = 100000

# prior theta
x = seq(0,1,length=1000)
theta_pri = dbeta(x,37.1,16.96)
# plot Jeffreys prior
plot(x,theta_pri,type="l",xlab = "theta",ylab="beta",main="Conjugate Prior ~ Beta(37.1,16.96)")
#### highest belief with 0/1, lowest with 0.5
# prior HPD
# hpd for theta
hpd_u1  = hdi(qbeta,shape1=37.1,shape2=16.96,credMass = 0.95)
#     lower     upper 
# 0.5627632 0.8059640  
#plot
# filters the values so that only the ones within the bounds are returned
bounds_filter1 <- x >= hpd_u1[1] & x <= hpd_u1[2]
x_within_bounds1 <- x[bounds_filter1]
theta_within_bounds1 <- theta_pri[bounds_filter1]

plot(x,theta_pri,type="l",xlab = "theta",ylab="beta",main="Prior 95% HPD for theta")
x_polygon1 <- c(hpd_u1[1], x_within_bounds1, hpd_u1[2])
theta_polygon1 <- c(0, theta_within_bounds1, 0)
polygon(x_polygon1, theta_polygon1, col = "skyblue")
text(0.5, 4, "0.56",  cex = 0.8, col = "gray30")
text(0.85, 4, "0.81 ", cex = 0.8, col = "gray30")



# posterior theta
theta_post = dbeta(x,173.1,100.96) #  
# posterior HPD
library(HDInterval)
# hpd for theta
hpd_u  = hdi(qbeta,shape1=173.1,shape2=100.96,credMass = 0.95)
#     lower     upper 
# 0.5744035 0.6882864 
#plot
# filters the values so that only the ones within the bounds are returned
bounds_filter <- x >= hpd_u[1] & x <= hpd_u[2]
x_within_bounds <- x[bounds_filter]
theta_within_bounds <- theta_post[bounds_filter]

plot(x,thate_post,type="l",xlab = "theta",ylab="beta",main="95% HPD for theta")
x_polygon <- c(hpd_u[1], x_within_bounds, hpd_u[2])
theta_polygon <- c(0, theta_within_bounds, 0)
polygon(x_polygon, theta_polygon, col = "skyblue")
text(0.5, 4, "0.574",  cex = 0.8, col = "gray30")
text(0.75, 4, "0.688 ", cex = 0.8, col = "gray30")


# preictive 
ypred = rbinom(n.sim,1,theta_post)
table(ypred)/n.sim
#       0       1 
# 0.38061 0.61939  
# plot predictive
hist(ypred,freq = F, col = "cornflowerblue",border = "wheat",right = FALSE, main="predictive value",xlab="",ylab=" ")
text(0, 8,"0.381",cex = 0.8, col = "black")
text(1, 12.6,"0.619",cex = 0.8, col = "black")
title(xlab = "")

# CDF(>0.5)
pbeta(0.5,173.1,100.96,lower.tail = F)
# 0.9999943

########################################## 3. compare 3 curves ########################################## 
source('~/Documents/Min/Study/Math/Statistics/SJSU MATH 264 Bayesian Analysis/project/264 project/beta_prior_post.R')
calcPosteriorForProportion(136, 220, 37.1,16.96)
# [1] "mode for prior= 0.693430656934307 , for likelihood= 0.618181818181818 , for posterior= 0.63258104829817"
# [1] "mean for prior= 0.686274509803922 , for likelihood= 0.617117117117117 , for posterior= 0.631613515288623"
# [1] "sd for prior= 0.0625324916629157 , for likelihood= 0.0325510004033116 , for posterior= 0.0290846493459548"

# check flat curves
calcPosteriorForProportion(136, 220, 1,1)
# [1] "mode for prior= NaN , for likelihood= 0.618181818181818 , for posterior= 0.618181818181818"
# [1] "mean for prior= 0.5 , for likelihood= 0.617117117117117 , for posterior= 0.617117117117117"
# [1] "sd for prior= 0.288675134594813 , for likelihood= 0.0325510004033116 , for posterior= 0.0325510004033116"


#-------------------
# model checking
#-------------------
################# sensitivity(dif conjugate)
calcPosteriorForProportion(136, 220, 35,100) # mode=0.26,var=0.0014
calcPosteriorForProportion(136, 220, 35,50) # mode=0.4,var=0.0028
calcPosteriorForProportion(136, 220, 35,35) # mode=0.5,var=0.0035
calcPosteriorForProportion(136, 220, 37.1,16.96) # mode=0.7,var=0.0039

################# flat/jef/conjugate
calcPosteriorForProportion(136, 220, 1,1) # flat
calcPosteriorForProportion(136, 220, 0.5,0.5) #jef
calcPosteriorForProportion(136, 220, 37.1,16.96) #conjugate

################# indenpendence
library(ggplot2)
library(tidyr)
set.seed(123)
# Compute test statistic for the data. Test statistic is the number of switches from 0 to 1 or from 1 to 0.
y = data$type
Ty <- sum(diff(y) != 0) + 0.0

# Sufficient statistics
n <- length(y)
s <- sum(y)

# Compute test statistic for the replicate data.
rb <- function(s, n) {
  p <- rbeta(1, s+37.1, n-s+16.96)
  yr <- rbinom(n, 1, p)
  sum(diff(yr) != 0) + 0.0
}
Tyr <- data.frame(x = replicate(n.sim, rb(s, n)))

# Compute posterior predictive p-value
mean(Tyr<=Ty)
# [1] 0.53846


# Plot test statistics for the data and replicates. Vertical line corresponds to the original data, and the histogram to the replicate data.
title1 <- 'Posterior Predictive Checking for Indenpendence
P-value = 0.538'
ggplot(data = Tyr) +
  geom_histogram(aes(x = x), fill = 'darkblue',color = 'black', binwidth = 1) +
  geom_vline(aes(xintercept = x), data = data.frame(x = Ty), color = 'red') +
  labs(x = '', y = '', title = title1) +
  scale_y_continuous(breaks=NULL)


################# mean
theta_rep = rbeta(n.sim,173.1,100.96)
yrep = mapply(rbinom,n=220,size=1,prob=theta_rep)
obs.mean <- mean(data$type)
sim.mean <- apply(yrep, 2, mean)
pval <- length(sim.mean[sim.mean >= obs.mean]) / n.sim
# 0.64398

hist(sim.mean,col = "cornflowerblue", border = "wheat",
     freq = FALSE, right = FALSE, breaks = 30,main="Posterior Predictive Checking for Mean",xlab=" ",ylab=" ")

lines(rep(obs.mean, 2), c(0, 10), col = "red", lwd = 1.2)

text(0.625, 10, "observed mean", srt = 90, cex = 0.8, col = "gray30")

title(xlab = "Replicated mean")
mtext(side = 3, paste0("p-value = ", round(pval, digits = 2)))



