#'Compute Power for Model 7 Conditional Processes Using Joint Significance
#'Requires correlations between all variables as sample size
#'Several values default to zero if no value provided
#'This is the recommended approach for determining power
#'@param rxy Correlation between DV (y) and predictor (x)
#'@param rxm Correlation between predictor (x) and mediator (m)
#'@param rxw Correlation between predictor (x) and moderator (w)
#'@param rxxw Correlation between predictor (x) and interaction term (xw) - defaults to 0
#'@param rmy Correlation between DV (y) and mediator (m)
#'@param rwy Correlation between DV (y) and moderator (w)
#'@param rwm Correlation between moderator (w) and mediator (m)
#'@param rxwy Correlation between DV (y) and interaction (xw) - defaults to 0
#'@param rwxw Correlation between moderator (w) and interaction (xw) - defaults to 0
#'@param rmxw Correlation between mediator (m) and interaction (xw) - Key value
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@param rep Number of samples drawn (defaults to 5000)
#'@examples \donttest{modmed7(rxm=.4, rxw=.2, rxy=.3, rwm=.2, rmxw=.1, rmy=.3,n=200)}
#'@return Power for Model 7 Conditional Processes
#'@export
#'
#'
modmed7<-function(rxm, rxw, rxxw=0, rxy,
                rwm, rwxw=0, rwy=0,
                rmxw, rmy, rxwy=0,
                alpha=.05,rep=1000,n=NULL){
V1<-NA;V2<-NA;V3<-NA;V4<-NA;V5<-NA

set.seed(1235)
out <- MASS::mvrnorm(100000, mu = c(0,0,0,0,0),
                     Sigma = matrix(c(1.0,rxw,rxm,rxxw,rxy,
                                      rxw,1.0,rwm,rwxw,rwy,
                                      rxm,rwm,1.0,rmxw,rmy,
                                      rxxw,rmxw,rwxw,1.0,rxwy,
                                      rxy,rwy,rmy,rxwy,1.0),
                                    ncol = 5),
                     empirical = TRUE)
out<-as.data.frame(out)

out<-dplyr::rename(out, x = V1, w= V2, m= V3, xw = V4, y = V5)


nruns = rep
a = numeric(nruns)
b = numeric(nruns)
sea = numeric(nruns)
seb = numeric(nruns)
pa3<-NA
pb<-NA


for (i in 1:nruns)
{
  samp <- out[ sample(nrow(out), n), ]
  test_a <- stats::lm(m ~ x + w + xw, data = samp)
  test_b <- stats::lm(y ~ x + m, data = samp)
  apath<-summary(test_a)
  bpath<-summary(test_b)
  pa3[i]<-apath$coefficients[4,4]
  pb[i]<-bpath$coefficients[3,4]
  power = data.frame(pa3=pa3,pb=pb)
  power$jsa["JointSiga"]<-NA
  power$jsa[pa3 < alpha] <- 1
  power$jsa[pa3 >= alpha] <- 0
  power$jsb1["JointSigb1"]<-NA
  power$jsb1[pb < alpha] <- 1
  power$jsb1[pb >= alpha] <- 0
  power$joint<-power$jsa*power$jsb
  JSa<-mean(power$joint)

}
message("Sample size is ",n)
message("Power for Index of Moderated Mediation ", JSa)
}


