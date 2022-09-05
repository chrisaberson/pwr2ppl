#'Compute Power for Conditional Process Model 14 Joint Significance
#'Requires correlations between all variables as sample size.
#'This is the recommended approach for determining power
#'@param rxy Correlation between DV (y) and predictor (x)
#'@param rxm Correlation between predictor (x) and mediator (m)
#'@param rxw Correlation between predictor (x) and moderator (w)
#'@param rxxw Correlation between predictor (x) and xweraction term (xw) - defaults to 0
#'@param rmy Correlation between DV (y) and mediator (m)
#'@param rwy Correlation between DV (y) and moderator (w)
#'@param rwm Correlation between moderator (w) and mediator (m)
#'@param rxwy Correlation between DV (y) and xweraction (xw) - defaults to 0
#'@param rxww Correlation between moderator (w) and xweraction (xw) - defaults to 0
#'@param rxwm Correlation between mediator (m) and xweraction (xw) - Key value
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@param rep Number of samples drawn (defaults to 5000)
#'@examples \donttest{modmed14(rxw=.2, rxm=.2, rxy=.31,rwy=.35, rxwy=.2,
#'rmy=.32, n=200, rep=1000,alpha=.05)}
#'@return Power for Model 14 Conditional Processes
#'@export
#'
#'

modmed14<-function(rxw, rxm, rxxw=0, rxy,
                  rwm=0, rxww=0, rwy,
                  rxwm=0, rxwy,
                  rmy, n,alpha=.05,rep=5000)
{
V1<-NA;V2<-NA;V3<-NA;V4<-NA;V5<-NA;V6<-NA
set.seed(1235)
out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0),
                     Sigma = matrix(c(1.0,rxw,rxm,rxxw,rxy,
                                      rxw,1.0,rwm,rxww,rwy,
                                      rxm,rwm,1.0,rxwm,rmy,
                                      rxxw,rxww,rxwm,1.0,rxwy,
                                      rxy,rwy,rmy,rxwy,1.0),
                                    ncol = 5),
                     empirical = TRUE)
out<-as.data.frame(out)

out<-dplyr::rename(out, x = V1, w= V2, m= V3, xw = V4, y = V5)


nruns = rep
a = numeric(nruns)
b = numeric(nruns)
pa1<-NA
pb3<-NA

for (i in 1:nruns)
{
  samp <- out[ sample(nrow(out), n), ]
  test_a <- stats::lm(m ~ x , data = samp)
  test_b <- stats::lm(y ~ x + m + w + xw, data = samp)
  apath<-summary(test_a)
  bpath<-summary(test_b)
  pa1[i]<-apath$coefficients[2,4] # x to m
  pb3[i]<-bpath$coefficients[5,4] #xw to y
  power = data.frame(pa1=pa1,pb3=pb3)

  power$jsa1["JoxwSiga1"]<-NA
  power$jsa1[pa1 < alpha] <- 1
  power$jsa1[pa1 >= alpha] <- 0
  power$jsb3["JoxwSigb1"]<-NA
  power$jsb3[pb3 < alpha] <- 1
  power$jsb3[pb3 >= alpha] <- 0
  power$joxw<-power$jsa1*power$jsb3
  JSa<-mean(power$joxw)
}


message("Sample size is ",n)
message("Power for Index of Moderated Mediation ", JSa)
}




