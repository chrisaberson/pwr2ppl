#'Compute Power for Serial Mediation Effects
#'Requires correlations between all variables as sample size.
#'This approach calculates power for the serial mediation using
#'joint significance (recommended)
#'@param rxy Correlation between DV (y) and predictor (x)
#'@param rxm1 Correlation between predictor (x) and first mediator (m1)
#'@param rxm2 Correlation between predictor (x) and second mediator (m2)
#'@param rym1 Correlation between DV (y) and first mediator (m1)
#'@param rym2 Correlation between DV (y) and second mediator (m2)
#'@param rm1m2 Correlation first mediator (m1) and second mediator (m2)
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@param reps Number of samples (default is 1000)
#'@examples
#'\donttest{medserial(rxm1=.3, rxm2=.3, rxy=-.35, rym1=-.5,rym2=-.5,
#'rm1m2=.7,n=150)}
#'@return Power for Serial Mediated (Indirect) Effects
#'@export
#'
#'

medserial<-function(rxm1,rxm2,rxy,rm1m2,rym1,rym2,n,alpha=.05, reps=1000)
{
  V1<-V2<-V3<-V4<-NA
  pop <- MASS::mvrnorm(100000, mu = c(0,0,0,0),
                       Sigma = matrix(c(1.0,rxm1,rxm2, rxy,
                                        rxm1,1.0,rm1m2, rym1,
                                        rxm2,rm1m2,1.0,rym2,
                                        rxy, rym1, rym2,1.0), ncol = 4),
                       empirical = TRUE)
  pop<-as.data.frame(pop)
  pop<-dplyr::rename(pop, x = V1, m1 = V2, m2 = V3, y = V4)
  set.seed(1234)
  nruns = reps
  a = numeric(nruns)
  b = numeric(nruns)
  pa1<-NA
  pa2<-NA
  pb1<-NA
  pb2<-NA
  pd<-NA
  for (i in 1:nruns)
  {
samp <- pop[ sample(nrow(pop), n), ]
pathsa1<-summary(lm(m1~x,data=samp))
pa1[i]<-pathsa1$coefficients[2,4]
pathsa2d<-summary(lm(m2~x+m1,data=samp))
pa2[i]<-pathsa2d$coefficients[2,4]
pd[i]<-pathsa2d$coefficients[3,4]
pathsb<-summary(lm(y~x+m1+m2,data=samp))
pb1[i]<-pathsb$coefficients[3,4]
pb2[i]<-pathsb$coefficients[4,4]


power = data.frame(pa1=pa1,pa2=pa2, pb1=pb1,pb2=pb2,pd=pd)
power$jsa1["JointSiga1"]<-NA
power$jsa1[pa1 < alpha] <- 1
power$jsa1[pa1 >= alpha] <- 0
power$jsa2["JointSiga2"]<-NA
power$jsa2[pa2 < alpha] <- 1
power$jsa2[pa2 >= alpha] <- 0
power$jsb1["JointSigb1"]<-NA
power$jsb1[pb1 < alpha] <- 1
power$jsb1[pb1 >= alpha] <- 0
power$jsb2["JointSigb2"]<-NA
power$jsb2[pb2 < alpha] <- 1
power$jsb2[pb2 >= alpha] <- 0
power$jsd["JointSigd"]<-NA
power$jsd[pd < alpha] <- 1
power$jsd[pd >= alpha] <- 0
power$jointm1<-power$jsa1*power$jsb1
power$jointm2<-power$jsa2*power$jsb2
power$jointm12<-power$jsa1*power$jsd*power$jsb2
JSm1<-mean(power$jointm1)
JSm2<-mean(power$jointm2)
JSm12<-mean(power$jointm12)}
message("Power for n = ", n,", mediator 1", " is ", JSm1)
message("Power for n = ", n,", mediator 2", " is ", JSm2)
message("Power for n = ", n,", serial mediation", " is ", JSm12)
result <- data.frame(matrix(ncol = 4))
colnames(result) <- c( "n","Power Mediator 1 (x-ma-mb)", "Power Mediator 2 (ma-mb-y)","Power Double (x-ma-mb-y)")
result[, 1]<-n
result[, 2]<-JSm1
result[, 3]<-JSm2
result[, 4]<-JSm12

output<-na.omit(result)
rownames(output)<- c()
invisible(output)

}
