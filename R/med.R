#'Compute Power for Mediated (Indirect) Effects
#'Requires correlations between all variables as sample size.
#'@param rxy Correlation between DV (y) and predictor (x)
#'@param rxm1 Correlation between predictor (x) and first mediator (m1)
#'@param rxm2 Correlation between predictor (x) and second mediator (m2)
#'@param rxm3 Correlation between predictor (x) and third mediator (m3)
#'@param rxm4 Correlation between predictor (x) and fourth mediator (m4)
#'@param rym1 Correlation between DV (y) and first mediator (m1)
#'@param rym2 Correlation between DV (y) and second mediator (m2)
#'@param rym3 Correlation DV (y) and third mediator (m3)
#'@param rym4 Correlation DV (y) and fourth mediator (m4)
#'@param rm1m2 Correlation first mediator (m1) and second mediator (m2)
#'@param rm1m3 Correlation first mediator (m1) and third mediator (m3)
#'@param rm1m4 Correlation first mediator (m1) and fourth mediator (m4)
#'@param rm2m3 Correlation second mediator (m2) and third mediator (m3)
#'@param rm2m4 Correlation second mediator (m2) and fourth mediator (m4)
#'@param rm3m4 Correlation third mediator (m3) and fourth mediator (m4)
#'@param n Sample size
#'@param mvars Number of Mediators
#'@param alpha Type I error (default is .05)
#'@examples
#'med(rxm1=.25, rxy=-.35, rym1=-.5,mvars=1, n=150)
#'med(rxm1=.3, rxm2=.3, rxm3=.25, rxy=-.35, rym1=-.5,rym2=-.5, rym3 = -.5,
#'rm1m2=.7, rm1m3=.4,rm2m3=.4, mvars=3, n=150)
#'@return Power for Mediated (Indirect) Effects
#'@export
#'
#'

med<-function(rxm1,rxm2=0, rxm3=0, rxm4=0,rxy,rym1, rym2=0,rym3=0, rym4=0, rm1m2=0,rm1m3=0,
                   rm1m4=0, rm2m3=0, rm2m4=0, rm3m4=0,alpha=.05,mvars,n)
{
  V1<-V2<-V3<-V4<-V5<-V6<-NA
  if(mvars==1){
    out <- MASS::mvrnorm(n, mu = c(0,0,0),
                   Sigma = matrix(c(1.0,rxm1,rxy,
                                    rxm1,1.0,rym1,
                                    rxy,rym1,1.0),
                                    ncol = 3),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, x = V1, m1 = V2, y = V3)
    set.seed(1234)
    model <- '
    # direct effect
    y ~ c*x
    # mediator
    m1 ~ a1*x
    y ~ b1*m1
    # indirect effect (a*b)
    ab1 := a1*b1
    # total effect (c-prime)
    total := c + (a1*b1)
    '
fitm <- lavaan::sem(model, data = out)
est<-lavaan::parameterestimates(fitm)
zab1<-abs(est$z[7])
tabled<-abs(stats::qnorm(alpha/2))
zpowerab1<-tabled-abs(zab1)
powerab1<-round(1-stats::pnorm(zpowerab1),4)
message("Power for n =", n,", mediator 1", " is ", powerab1)
result <- data.frame(matrix(ncol = 2))
colnames(result) <- c( "n","Power Mediator 1")
result[, 1]<-n
result[, 2]<-powerab1
output<-na.omit(result)
rownames(output)<- c()
}
if(mvars==2){
  out <- MASS::mvrnorm(n, mu = c(0,0,0,0),
               Sigma = matrix(c(1.0,rxm1,rxm2, rxy,
                                rxm1,1.0,rm1m2, rym1,
                                rxm2,rm1m2,1.0,rym2,
                                rxy, rym1, rym2,1.0), ncol = 4),
                                empirical = TRUE)
out<-as.data.frame(out)
out<-dplyr::rename(out, x = V1, m1 = V2, m2 = V3, y = V4)
set.seed(1234)
model <- '
# direct effect
y ~ c*x
# mediator
m1 ~ a1*x
m2 ~ a2*x
y ~ b1*m1
y ~ b2*m2
# indirect effect (a*b)
ab1 := a1*b1
ab2 := a2*b2
totalind := ab1+ab2
# total effect (c-prime)
total := c + (a1*b1)+(a2*b2)
#Cors
m1~~m2
'
fitm <- lavaan::sem(model, data = out)
est<-lavaan::parameterestimates(fitm)
zab1<-abs(est$z[11])
zab2<-abs(est$z[12])
zall<-abs(est$z[13])
tabled<-abs(stats::qnorm(alpha/2))
zpowerab1<-tabled-abs(zab1)
zpowerab2<-tabled-abs(zab2)
zpowerall<-tabled-abs(zall)
powerab1<-round(1-stats::pnorm(zpowerab1),4)
powerab2<-round(1-stats::pnorm(zpowerab2),4)
powerall<-round(1-stats::pnorm(zpowerall),4)
message("Power for n = ", n,", mediator 1", " is ", powerab1)
message("Power for n = ", n,", mediator 2", " is ", powerab2)
message("Power for n = ",n,", Total Mediation", " is ", powerall)
result <- data.frame(matrix(ncol = 4))
colnames(result) <- c( "n","Power Mediator 1", "Power Mediator 2", "Power Total")
result[, 1]<-n
result[, 2]<-powerab1
result[, 3]<-powerab2
result[, 4]<-powerall

output<-na.omit(result)
rownames(output)<- c()
}

else if(mvars==3){
  out <- MASS::mvrnorm(n, mu = c(0,0,0,0,0),
                 Sigma = matrix(c(1.0,rxm1,rxm2, rxm3, rxy,
                                  rxm1,1.0,rm1m2,rm1m3,rym1,
                                  rxm2,rm1m2,1.0,rm2m3,rym2,
                                  rxm3,rm1m3,rm2m3,1.0,rym3,
                                  rxy, rym1, rym2,rym3,1.0), ncol = 5),
                 empirical = TRUE)
  out<-as.data.frame(out)
  out<-dplyr::rename(out, x = V1, m1 = V2, m2 = V3, m3=V4, y = V5)
  set.seed(1234)
  model <- '
  # direct effect
  y ~ c*x
  # mediator
  m1 ~ a1*x
  m2 ~ a2*x
  m3 ~ a3*x
  y ~ b1*m1
  y ~ b2*m2
  y ~ b3*m3
  # indirect effect (a*b)
  ab1 := a1*b1
  ab2 := a2*b2
  ab3 := a3*b3
  totalind := ab1+ab2+ab3
  # total effect (c-prime)
  total := c + (a1*b1)+(a2*b2)+(a3*b3)
  #Cors
  m1~~m2
  m1~~m3
  m2~~m3
  '
  fitm <- lavaan::sem(model, data = out)
  est<-lavaan::parameterestimates(fitm)
  zab1<-abs(est$z[16])
  zab2<-abs(est$z[17])
  zab3<-abs(est$z[18])
  zall<-abs(est$z[19])
  tabled<-abs(stats::qnorm(alpha/2))
  zpowerab1<-tabled-abs(zab1)
  zpowerab2<-tabled-abs(zab2)
  zpowerab3<-tabled-abs(zab3)
  zpowerall<-tabled-abs(zall)
  powerab1<-round(1-stats::pnorm(zpowerab1),4)
  powerab2<-round(1-stats::pnorm(zpowerab2),4)
  powerab3<-round(1-stats::pnorm(zpowerab3),4)
  powerall<-round(1-stats::pnorm(zpowerall),4)
  message("Power for n = ", n,", mediator 1", " is ", powerab1)
  message("Power for n = ", n,", mediator 2", " is ", powerab2)
  message("Power for n =", n,", mediator 3", " is ", powerab3)
  message("Power for n = ",n,", Total Mediation", " is ", powerall)
  result <- data.frame(matrix(ncol = 5))
  colnames(result) <- c( "n","Power Mediator 1", "Power Mediator 2", "Power Mediator 3", "Power Total")
  result[, 1]<-n
  result[, 2]<-powerab1
  result[, 3]<-powerab2
  result[, 4]<-powerab3
  result[, 5]<-powerall
  output<-na.omit(result)
  rownames(output)<- c()
  est
}

  else if(mvars==4){
    out <- MASS::mvrnorm(n, mu = c(0,0,0,0,0,0),
                   Sigma = matrix(c(1.0,rxm1,rxm2, rxm3, rxm4, rxy,
                                    rxm1,1.0,rm1m2,rm1m3,rm1m4,rym1,
                                    rxm2,rm1m2,1.0,rm2m3,rm2m4,rym2,
                                    rxm3,rm1m3,rm2m3,1.0,rm3m4,rym3,
                                    rxm4,rm1m4,rm2m4,rm3m4,1.0,rym4,
                                    rxy, rym1, rym2,rym3,rym4,1.0), ncol = 6),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, x = V1, m1 = V2, m2 = V3, m3=V4, m4=V5,y = V6)
    set.seed(1234)
    model <- '
    # direct effect
    y ~ c*x
    # mediator
    m1 ~ a1*x
    m2 ~ a2*x
    m3 ~ a3*x
    m4 ~ a4*x
    y ~ b1*m1
    y ~ b2*m2
    y ~ b3*m3
    y ~ b4*m4
    # indirect effect (a*b)
    ab1 := a1*b1
    ab2 := a2*b2
    ab3 := a3*b3
    ab4 := a4*b4
    totalind := ab1+ab2+ab3+ab4
    # total effect (c-prime)
    total := c + (a1*b1)+(a2*b2)+(a3*b3)+(a4*b4)
    #Cors
    m1~~m2
    m1~~m3
    m1~~m4
    m2~~m3
    m2~~m4
    m3~~m4
    '
    fitm <- lavaan::sem(model, data = out)
    est<-lavaan::parameterestimates(fitm)
    zab1<-est$z[22]
    zab2<-est$z[23]
    zab3<-est$z[24]
    zab4<-est$z[25]
    zall<-est$z[26]
    tabled<-abs(stats::qnorm(alpha/2))
    zpowerab1<-tabled-abs(zab1)
    zpowerab2<-tabled-abs(zab2)
    zpowerab3<-tabled-abs(zab3)
    zpowerab4<-tabled-abs(zab4)
    zpowerall<-tabled-abs(zall)
    powerab1<-round(1-stats::pnorm(zpowerab1),4)
    powerab2<-round(1-stats::pnorm(zpowerab2),4)
    powerab3<-round(1-stats::pnorm(zpowerab3),4)
    powerab4<-round(1-stats::pnorm(zpowerab4),4)
    powerall<-round(1-stats::pnorm(zpowerall),4)
    message("Power for n = ", n,", mediator 1", " is ", powerab1)
    message("Power for n = ", n,", mediator 2", " is ", powerab2)
    message("Power for n = ", n,", mediator 3", " is ", powerab3)
    message("Power for n = ", n,", mediator 4", " is ", powerab4)
    message("Power for n = ",n,", Total Mediation", " is ", powerall)
    result <- data.frame(matrix(ncol = 6))
    colnames(result) <- c( "n","Power Mediator 1", "Power Mediator 2", "Power Mediator 3","Power Mediator 4", "Power Total")
    result[, 1]<-n
    result[, 2]<-powerab1
    result[, 3]<-powerab2
    result[, 4]<-powerab3
    result[, 5]<-powerab4
    result[, 6]<-powerall
    output<-na.omit(result)
    rownames(output)<- c()
    }
invisible(output)

}
