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
#  0   1 
# 84 136 
table(data$type)/220
#         0         1 
# 0.3818182 0.6181818 



#-------------------
# predictive
# 
# Jef Prior: Beta(0.5,0.5)
# sampling: 
# posterior:beta(0.5+sum(y), n-sum(y)+0.5) = beta(136.5,84.5); mode=135.5/219=0.6187215
# predict:bin(1,1,theta)

# #-------------------

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
y_pri = dbeta(x,0.5,0.5)
# plot Jeffreys prior
plot(x,y_pri,type="l",xlab = "theta",ylab="beta",main="Jeffreys Prior Distribution for theta")
#### highest belief with 0/1, lowest with 0.5


# posterior theta
y_post = dbeta(x,136.5,84.5) #  
# HPD
library(HDInterval)
# hpd for theta
hpd_u  = hdi(qbeta,shape1=136.5,shape2=84.5,credMass = 0.95)
#     lower     upper 
# 0.5535067 0.6811925 
#plot
# filters the values so that only the ones within the bounds are returned
bounds_filter <- x >= hpd_u[1] & x <= hpd_u[2]
x_within_bounds <- x[bounds_filter]
y_within_bounds <- y_post[bounds_filter]

plot(x,y_post,type="l",xlab = "theta",ylab="beta",main="Posterior Distribution for theta")
x_polygon <- c(hpd_u[1], x_within_bounds, hpd_u[2])
y_polygon <- c(0, y_within_bounds, 0)
polygon(x_polygon, y_polygon, col = "skyblue")
text(hpd_u[1], 4, "0.554", srt = 90, cex = 0.8, col = "gray30")
text(hpd_u[2], 4, "0.681 ", srt = 90, cex = 0.8, col = "gray30")


# preictive 
ypred = rbinom(n.sim,1,theta_post)
table(ypred)/n.sim
#       0       1 
# 0.38179 0.61821 
# plot predictive
hist(ypred,freq = F, col = "cornflowerblue",border = "wheat",right = FALSE, main="predictive value",xlab="",ylab=" ")
text(0, 8,"0.381",cex = 0.8, col = "black")
text(1, 12.6,"0.618",cex = 0.8, col = "black")
title(xlab = "")


########################################## 3. compare 3 curves ########################################## 
source('~/Documents/Min/Study/Math/Statistics/SJSU MATH 264 Bayesian Analysis/project/264 project/beta_prior_post.R')
calcPosteriorForProportion(136, 220, 0.5, 0.5)
# [1] "mode for prior= 0.5 , for likelihood= 0.618181818181818 , for posterior= 0.618721461187215"
# [1] "mean for prior= 0.5 , for likelihood= 0.617117117117117 , for posterior= 0.617647058823529"
# [1] "sd for prior= 0.353553390593274 , for likelihood= 0.0325510004033116 , for posterior= 0.0326156410793778"
#-------------------
# model checking
#-------------------
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
  p <- rbeta(1, s+0.5, n-s+0.5)
  yr <- rbinom(n, 1, p)
  sum(diff(yr) != 0) + 0.0
}
Tyr <- data.frame(x = replicate(n.sim, rb(s, n)))

# Compute posterior predictive p-value
mean(Tyr<=Ty)
# [1] 0.47268


# Plot test statistics for the data and replicates. Vertical line corresponds to the original data, and the histogram to the replicate data.
title1 <- 'Binomial example - number of changes?
Pr(T(yrep,theta) <= T(y,theta)|y) = 0.473'
ggplot(data = Tyr) +
  geom_histogram(aes(x = x), fill = 'darkblue',color = 'black', binwidth = 1) +
  geom_vline(aes(xintercept = x), data = data.frame(x = Ty), color = 'red') +
  labs(x = '', y = '', title = title1) +
  scale_y_continuous(breaks=NULL)




