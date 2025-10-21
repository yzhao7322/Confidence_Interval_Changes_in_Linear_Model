library('tseries')
library('zoo')
library('pracma')
library('CPAT')
library('cointReg')
library('sde')
library('MASS')
library('stats')
library('pracma')
library('fGarch')
library('stepR')
library('strucchange')

# load confidence_interval_functions.R
source("confidence_interval_functions.R")
# load ilr_code.R
source("ilr_code.R")

##################
### Application 1 - Section 6.
##################
# Slack and Cyclically Sensitive Inflation, JMCB, Stock and Watson, 2021.
# New Keynesian Phillips curve
pcurve = read.csv("your_working_path\\application1\\pcurve.csv",header=TRUE)
monthly_data = pcurve[,2:5] 
T_month = dim(monthly_data)[1] 


####################
# pre-analysis
####################
# y & x
infla_headline_month = infla_core_month = c(rep(NA, T_month-12) )
for (i in 13:T_month){
  # Year-over-year change in the rate of inflation 
	infla_headline_month[i-12] = (monthly_data[i,1]-monthly_data[i-12,1])/monthly_data[i-12,1]*100
	infla_core_month[i-12] = (monthly_data[i,2]-monthly_data[i-12,2])/monthly_data[i-12,2]*100
}
  # the four-quarter average of the CBO unemployment gap.
unemploy_month = rollmean(c(monthly_data[,4]-monthly_data[,3]),k=13)
T =  length(infla_core_month)


adf.test(infla_headline_month)
adf.test(infla_core_month)
adf.test(unemploy_month)
 

# core regression
rcore_month = lm(infla_core_month ~ unemploy_month -1 )
summary(rcore_month)
rcore_month_resid = residuals(rcore_month)

# Coefficients:
#                Estimate Std. Error t value Pr(>|t|)
# unemploy_month -0.61474    0.09168  -6.705 4.01e-11 ***


###### HeteroVar error check
wl=59
resv = ae = c(seq(1:(T-wl)))
for (i in 1:(T-wl)){
  fit <- lm(infla_core_month[i:(i+wl)] ~ unemploy_month[i:(i+wl)] -1) 
  resv[i] = sd(residuals(fit))
}
fit1 <- lm(infla_core_month ~ unemploy_month -1) 
ae = abs(residuals(fit1)[(T-(T-wl)+1):T])

resv1 =c(rep(NA,60),resv)

# plot variables and residuals
par(mfrow=c(1,2),mar=c(3,3,3,3))
plot(infla_core_month,type="l",cex=1.8,col='blue',lwd=3,xaxt="n",ylim=c(-5,15),ylab="Value",xlab="Date",main="Variables",cex.main=2.5,cex.axis = 2.5)
axis(1,at=seq(1,737,184),cex.axis=2.2,labels = c("Jan-1961", "May-1976", "Sep-1991", "Jan-2007", "May-2022"))
lines(unemploy_month,type="l",lwd=3,cex=1.8,col='red')
abline(h=0, col="springgreen4",lwd=2.5,lty=3)
legend("topright",legend=c("Inflation","Unemployment"),lty = c(1,1),lwd=3, col=c("blue","red"),cex=1.8,y.intersp=1.8,bty = "n",seg.len=0.5)
plot(resv,type="l",cex=1.8,col='black',lwd=3,xaxt="n",ylim=c(-5,15),ylab="Value",xlab="Date",main="",cex.main=2.5,cex.axis = 2.5)
# lines(ae,lwd=3,col='darkgoldenrod3')
axis(1,at=c(seq(1,678,169)[1:4],678),cex.axis=2.2,labels = c("Jan-1966", "Feb-1980", "Mar-1994", "Apr-2008", "May-2022"))
abline(h=0, col="springgreen4",lwd=2.5,lty=3)
# legend("topleft",legend=c("Absolute errors","SD of rolling window errors"),lty = c(1,1),lwd=3, col=c("black","darkgoldenrod3"),cex=1.8,y.intersp=1.8,bty = "n",seg.len=0.5)
legend("topleft",legend=c("Rolling window estimate of the SD of the residuals"),lty = c(1),lwd=3, col=c("black"),cex=1.8,y.intersp=1.8,bty = "n",seg.len=0.5)




###############################
# estimating changes and confidence intervals
###############################

### K_N^{SMO} method
y = as.matrix(infla_core_month)
x = as.matrix(unemploy_month)

# manually


# auto
bre_out = greedy_breaks(y, x, kappa=1/2, max_breaks = 7, min_seg_n = 20)
bre_out$breaks
res_breaks <- combinations_with_anchor(bre_out$breaks)
breaks_out = break_selection(y,x,res_breaks)
kappa = 1/2 
N = length(y)
### change point refinement
refined_br = change_point_refinement(y,x,breaks_out, kappa)
refined_br
# 170 260 330 427 
# 1975/02, 1982/08, 1988/06,  1996/07

###################
# 711 [2020/03] the break around Covid-19 at was not selected via BIC criteria
###################

# confidence interval
set.seed(55)
ksc_confidence_interval(y, x, refined_br, band_width=floor(sqrt(N)), kappa = 1/2, alpha = 0.95)
# [[1]]
# 97.5%  2.5% 
#   162   198
# 1974/6 - 1977/6
# [[2]]
# 97.5%  2.5%
#   254   283
# 1982/02 - 1984/7
# [[3]]
# 97.5%  2.5%
#   307   381
# 1986/07 - 1992/9
# [[4]]
# 97.5%  2.5%
#   424   428
# 1996/04 - 1996/08

### BP method
alpha = 0.95
y = as.matrix(infla_core_month)
x = as.matrix(unemploy_month)
bp.yx <- breakpoints(y ~ x)
bp.yx
#breaks  157 268 378 497
#    1974/01, 1983/04, 1992/06, 2002/05
confint(bp.yx,level=alpha)$confint 
#   2.5 % breakpoints 97.5 %
# 1   156 (Dec1973)       157    158 (Feb1974)
# 2   267 (Mar1983)       268    270 (Jun1983)
# 3   377 (May1992)       378    380 (Aug1992)
# 4   488 (Aug2001)       497    499 (Jul2002)


# ILR method
y = as.matrix(infla_core_month)
x = as.matrix(unemploy_month)

ilr.yx <- ilr2015(y, x)
ilr.yx
#breaks   157  268  378   497
# lower - 156  266  366   447
# upper - 159  277  382   549 (700), noticed that 700 here is irrelavant, if checking "results$csm95" in ilr.yx function, it explains.
# (Dec1973) (Mar1974)
# (Feb1983) (Jan1984)
# (Jun1991) (Oct1992)
# (Mar1998) (Sep2006)


## plot of change point + confidence interval
#  K_N^{SMO}
col1="gray70"
plot(infla_core_month,type="l",cex=1.8,col='blue',lwd=3,xaxt="n",ylim=c(-5,15),ylab="",xlab="",main="",cex.main=2.5,cex.axis = 2.5)
rect(162, -10, 198, 20, border = col1, col = col1)
rect(254, -10, 283, 20, border = col1, col = col1)
rect(307, -10, 381, 20, border = col1, col = col1)
rect(424, -10, 428, 20, border = col1, col = col1)
lines(infla_core_month,type="l",cex=1.8,col='blue',lwd=3)
axis(1,at=seq(1,737,184),cex.axis=2.2,labels = c("Jan-1961", "May-1976", "Sep-1991", "Jan-2007", "May-2022"))
lines(unemploy_month,type="l",lwd=3,cex=1.8,col='red')
lines(resv1,lty=2,cex=2,col='magenta1',lwd=4)
abline(h=0, col="springgreen4",lwd=2.5,lty=3)
abline(v=c(170,260,330,427), col="black",lwd=4,lty=1)
legend(x = c(480, 480), y = c(15, 15),legend=c("Inflation","Unemployment","Rolling window SD","Breaks"),lty = c(1,1,1,1),lwd=3, col=c("blue","red","magenta1","black"),cex=1.8,y.intersp=1.2,bty = "n",seg.len=0.5)
 
# BP
col1="gray70"
plot(infla_core_month,type="l",cex=1.8,col='blue',lwd=3,xaxt="n",ylim=c(-5,15),ylab="",xlab="",main="",cex.main=2.5,cex.axis = 2.5)
rect(156, -10, 158, 20, border = col1, col = col1)
rect(267, -10, 270, 20, border = col1, col = col1)
rect(377, -10, 380, 20, border = col1, col = col1)
rect(488, -10, 499, 20, border = col1, col = col1)
lines(infla_core_month,type="l",cex=1.8,col='blue',lwd=3)
axis(1,at=seq(1,737,184),cex.axis=2.2,labels = c("Jan-1961", "May-1976", "Sep-1991", "Jan-2007", "May-2022"))
lines(unemploy_month,type="l",lwd=3,cex=1.8,col='red')
lines(resv1,lty=2,cex=2,col='magenta1',lwd=4)
abline(h=0, col="springgreen4",lwd=2.5,lty=3)
abline(v=c(157,268,378,497), col="black",lwd=4,lty=1)
legend(x = c(480, 480), y = c(15, 15),legend=c("Inflation","Unemployment","Rolling window SD","Breaks"),lty = c(1,1,1,1),lwd=3, col=c("blue","red","magenta1","black"),cex=1.8,y.intersp=1.2,bty = "n",seg.len=0.5)
 
# ilr2015
col1="gray70"
plot(infla_core_month,type="l",cex=1.8,col='blue',lwd=3,xaxt="n",ylim=c(-5,15),ylab="",xlab="",main="",cex.main=2.5,cex.axis = 2.5)
rect(156, -10, 159, 20, border = col1, col = col1)
rect(266, -10, 277, 20, border = col1, col = col1)
rect(366, -10, 382, 20, border = col1, col = col1)
rect(447, -10, 549, 20, border = col1, col = col1) 
lines(infla_core_month,type="l",cex=1.8,col='blue',lwd=3)
axis(1,at=seq(1,737,184),cex.axis=2.2,labels = c("Jan-1961", "May-1976", "Sep-1991", "Jan-2007", "May-2022"))
lines(unemploy_month,type="l",lwd=3,cex=1.8,col='red')
lines(resv1,lty=2,cex=2,col='magenta1',lwd=4)
abline(h=0, col="springgreen4",lwd=2.5,lty=3)
abline(v=c(157, 268, 378, 497), col="black",lwd=4,lty=1)
legend(x = c(480, 480), y = c(15, 15),legend=c("Inflation","Unemployment","Rolling window SD","Breaks"),lty = c(1,1,1,1),lwd=3, col=c("blue","red","magenta1","black"),cex=1.8,y.intersp=1.2,bty = "n",seg.len=0.5)



 
## estimate the slope of the Phillips curve
# Table 6.1

# BP / ILR
y1 = infla_core_month[1:157]
x1 = unemploy_month[1:157]
y2 = infla_core_month[158:268]
x2 = unemploy_month[158:268]
y3 = infla_core_month[269:378]
x3 = unemploy_month[269:378]
y4 = infla_core_month[379:497]
x4 = unemploy_month[379:497]
y5 = infla_core_month[498:737]
x5 = unemploy_month[498:737]

table2 = matrix(NA, 5, 2)
table2[1,1] = summary(lm(y1 ~ x1 -1 ))[4]$coefficients[1]
table2[1,2] = sd(residuals(lm(y1 ~ x1 -1 )))
table2[2,1] = summary(lm(y2 ~ x2 -1 ))[4]$coefficients[1]
table2[2,2] = sd(residuals(lm(y2 ~ x2 -1 )))
table2[3,1] = summary(lm(y3 ~ x3 -1 ))[4]$coefficients[1]
table2[3,2] = sd(residuals(lm(y3 ~ x3 -1 )))
table2[4,1] = summary(lm(y4 ~ x4 -1 ))[4]$coefficients[1]
table2[4,2] = sd(residuals(lm(y4 ~ x4 -1 )))
table2[5,1] = summary(lm(y5 ~ x5 -1 ))[4]$coefficients[1]
table2[5,2] = sd(residuals(lm(y5 ~ x5 -1 )))
round(table2,2)


# K_N^{SMO}
y1 = infla_core_month[1:170]
x1 = unemploy_month[1:170] 
y2 = infla_core_month[171:260]
x2 = unemploy_month[171:260]
y3 = infla_core_month[261:330]
x3 = unemploy_month[261:330] 
y4 = infla_core_month[331:427]
x4 = unemploy_month[331:427] 
y5 = infla_core_month[428:737]
x5 = unemploy_month[428:737]

table2 = matrix(NA, 5, 2)
table2[1,1] = summary(lm(y1 ~ x1 -1 ))[4]$coefficients[1]
table2[1,2] = sd(residuals(lm(y1 ~ x1 -1 )))
table2[2,1] = summary(lm(y2 ~ x2 -1 ))[4]$coefficients[1]
table2[2,2] = sd(residuals(lm(y2 ~ x2 -1 )))
table2[3,1] = summary(lm(y3 ~ x3 -1 ))[4]$coefficients[1]
table2[3,2] = sd(residuals(lm(y3 ~ x3 -1 )))
table2[4,1] = summary(lm(y4 ~ x4 -1 ))[4]$coefficients[1]
table2[4,2] = sd(residuals(lm(y4 ~ x4 -1 )))
table2[5,1] = summary(lm(y5 ~ x5 -1 ))[4]$coefficients[1]
table2[5,2] = sd(residuals(lm(y5 ~ x5 -1 )))
round(table2,2)


# post-Covid19

y6 = infla_core_month[711:737]
x6 = unemploy_month[711:737]
summary(lm(y6 ~ x6 -1 ))[4]$coefficients[1]
sd(residuals(lm(y6 ~ x6 -1 )))


##################
### Application 2
##################
########## Cryptocurrency
crypdat = read.csv("your_working_path\\application2\\crypfactors.csv",header=TRUE)
crypdat = read.csv("E:\\Working projects\\TimeofChange\\tex\\ET_sub\\RR\\code\\github\\crypfactors.csv",header=TRUE)
xfactor = as.matrix(crypdat[,2:4])
ydata = as.matrix(crypdat[,5])
rfdata = as.matrix(crypdat[,6])
yret = ydata-rfdata
T = length(yret) 

y = as.matrix(yret)
N = length(y)
x = as.matrix(cbind(rep(1,N),xfactor))
bre_out = greedy_breaks(y, x, kappa=1/2, max_breaks = 2, min_seg_n = 20)
res_breaks <- combinations_with_anchor(bre_out$breaks)
breaks_out = break_selection(y,x,res_breaks)
kappa = 1/2 

### change point refinement
refined_br = change_point_refinement(y,x,breaks_out, kappa)
refined_br
# 196 316 
# 2017/week 43 , 2020/week 07

# k_N^{(SMO)}(1/2)
ksc_confidence_interval(y, x, refined_br , band_width=floor(sqrt(N)), kappa = 1/2, alpha = 0.95)
# [[1]]
# 97.5%  2.5% 
#   196   197
# [[2]]
# 97.5%  2.5%
#   316   325

# Figure E.1

col1="gray70"
x<-1:339
plot(xfactor[,1],type="l",cex=1.8,col='red',lwd=2.5,xaxt ="n", yaxt ="n", ylim=c(-1.2,1.2),ylab="",xlab="",cex.axis = 2)
axis(side=4)
rect(196, -10, 197, 20, border = col1, col = col1)
rect(316, -10, 325, 20, border = col1, col = col1)
lines(xfactor[,1],type="l",cex=1.8,col='red',lwd=2.5)
lines(xfactor[,2],type="l",cex=1.8,col='cyan1',lwd=2.5)
lines(xfactor[,3],type="l",cex=1.8,col='goldenrod1',lwd=2.5)
axis(1,at=c(seq(1,339,84)[1:4],339),cex.axis=1.8,labels = c("20-Jan-2014", "24-Aug-2015", "10-Apr-2017", "19-Nov-2018", "20-Jul-2020"))
par(new=TRUE)
plot(x,yret,type="l",cex=1.8,col='blue',lwd=4,xaxt="n",ylim=c(-5,5),ylab="Value",xlab="Date",main="",cex.main=2.5,cex.axis = 2)
legend(x = c(2, 2), y = c(-0.5, -0.5),legend=c("Bitcoin","Cmkt","Csize","Cmom","Break"),lty = c(1,1,1,1,1),lwd=3, col=c("blue","red","cyan1","goldenrod1","black"),cex=1.8,y.intersp=1.2,bty = "n",seg.len=0.5)
abline(h=0, col="springgreen4",lwd=2.5,lty=3)
abline(v=c(196), col="black",lwd=4,lty=1)
abline(v=c(316), col="black",lwd=4,lty=1)


# results in Table E.1
regcryp1 = lm(yret[1:196] ~ xfactor[1:196,])
summary(regcryp1)[4]$coefficients
sd(residuals(regcryp1))

regcryp2 = lm(yret[197:316] ~ xfactor[197:316,])
summary(regcryp2)[4]$coefficients
sd(residuals(regcryp2))

regcryp3 = lm(yret[317:T] ~ xfactor[317:T,])
summary(regcryp3)[4]$coefficients
sd(residuals(regcryp3))
# estimating changes and confidence intervals
