#'Compute Power for Serial Mediation Effects
#'Requires correlations between all variables as sample size.
#'This approach calculates power for the serial mediation using
#'joint significance (recommended) and path coefficients
#'@param a1 path between predictor and first mediator
#'@param a2 path between predictor and first mediator
#'@param b1 Path between first mediator and dependent variable
#'@param b2 Path between first mediator and dependent variable
#'@param cprime Path between predictor and dependent variable
#'@param d Path first mediator (m1) and second mediator (m2)
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@param reps number of repetitions (1000 is default)
#'@examples
#'\donttest{medserial_paths(a1=.3, a2=.3, b1=.35,
#'b2=.3,d=.2,cprime=.1,n=150)}
#'@return Power for Serial Mediated (Indirect) Effects
#'@export
#'
#'

medserial_paths<-function(a1,a2,b1,b2,d,cprime,n,alpha=.05, reps=1000)
{
  V1<-NA;V2<-NA;V3<-NA;V4<-NA
  rxm1 <- a1 #x to M1
  rxm2 <- a2 + d*a1 #x to m2
  rm1m2 <- d + a1*a2 #m1 -> m2
  rxy <- cprime + a1*b1 + rxm2 #x to y
  rym1 <- a1*cprime + b1 + b2*rxm2 #m1 to y
  rym2 <- a2*rxm2 + b2 + b1*rm1m2 #m2 to y



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
