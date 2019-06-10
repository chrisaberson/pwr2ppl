#'Power for Comparing Dependent Coefficients in Multiple Regression with Two or Three Predictors
#'Requires correlations between all variables as sample size. Means, sds, and alpha are option. Also computes Power(All)
#'@param ry1 Correlation between DV (y) and first predictor (1)
#'@param ry2 Correlation between DV (y) and second predictor (2)
#'@param ry3 Correlation between DV (y) and third predictor (3)
#'@param r12 Correlation between first (1) and second predictor (2)
#'@param r13 Correlation between first (1) and third predictor (3)
#'@param r23 Correlation between second (2) and third predictor (3)
#'@param n Sample size
#'@param alpha Type I error (default is .05)
#'@examples
#'depb(ry1=.40, ry2=.40, ry3=-.40, r12=-.15, r13=-.60, r23=.25,n=110, alpha=.05)
#'@return Power for Comparing Dependent Coefficients in Multiple Regression with Two or Three Predictors
#'@export
#'
#'

depb<-function(ry1, ry2, ry3=NULL, r12, r13=NULL, r23=NULL,n=NULL, alpha=.05)
{
  pred<-NA
  pred[is.null(r23)]<-2
  pred[!is.null(r23)]<-3


  if (pred=="2")
  {pop <- MASS::mvrnorm(n, mu<-c(0,0,0), Sigma<-matrix(c(1, ry1, ry2,
                                                    ry1, 1, r12,
                                                    ry2, r12, 1),
                                                  ncol=3), empirical=TRUE)
  pop1<-data.frame(pop)
  values<-stats::lm(X1~X2+X3, pop1)
  values<-summary(values)
  b1<-(values$coefficients)[2,1] #grabs b from each analysis
  b2<-(values$coefficients)[3,1]
  seb1<-(values$coefficients)[2,2]
  seb2<-(values$coefficients)[3,2]
  df<-n-pred
  mat<-cbind(c(1,r12),c(r12,1))
  inv<-solve(mat)*mat
  pij<-inv[1,2] #inv of cor between pred of interest
  pii<-inv[1,1] #inv of cov, v1
  pjj<-inv[2,2] #inv of cov, v2
  den1<-seb1^2+seb2^2
  den2<-2*seb1*seb2
  den3<-pij/(pii+pjj)
  den<-(den1-(den2*den3))^.5
  t<-abs((abs(b1)-(abs(b2)))) / den
  lambda<-t^2
  df<-n-3
  minusalpha<-1-alpha
  Fb<-stats::qf(minusalpha, 1, df)
  power1<-round(1-stats::pf(Fb, 1,df,lambda),3)
  message("Sample size is ",n)
  message("Power Comparing b1 and b2 = ", power1)
  result <- data.frame(matrix(ncol = 2))
  colnames(result) <- c( "n","Power")
  result[, 1]<-n
  result[, 2]<-power1
  output<-na.omit(result)
  rownames(output)<- c()
    }

  if (pred=="3")
  {
    pop <- MASS::mvrnorm(n, mu<-c(0, 0, 0, 0), Sigma<-matrix(c(1, ry1, ry2, ry3,
                                                         ry1, 1, r12, r13,
                                                         ry2, r12,1, r23,
                                                         ry3, r13, r23, 1),
                                                         ncol=4), empirical=TRUE)

pop1<-data.frame(pop)
values<-stats::lm(X1~X2+X3+X4, pop1)
values<-summary(values)
b1<-(values$coefficients)[2,1] #grabs b from each analysis
b2<-(values$coefficients)[3,1]
b3<-(values$coefficients)[4,1]
seb1<-(values$coefficients)[2,2]
seb2<-(values$coefficients)[3,2]
seb3<-(values$coefficients)[4,2]
df<-n-pred
mat<-cbind(c(1,r12,r13),c(r12,1, r23),c(r13,r23,1))
inv<-solve(mat)*mat
# 1 vs 2
pij1<-inv[1,2] #inv of cor between pred of interest 1 vs. 2
pii1<-inv[1,1] #inv of cov, v1
pjj1<-inv[2,2] #inv of cov, v2
den1a<-seb1^2+seb2^2
den2a<-2*seb1*seb2
den3a<-pij1/(pii1+pjj1)
dena<-(den1a-(den2a*den3a))^.5
ta<-abs(abs(b1)-abs(b2))/ dena
lambdaa<-ta^2
#1 vs 3
pij2<-inv[1,3] #inv of cor between pred of interest 1 vs. 2
pii2<-inv[1,1] #inv of cov, v1
pjj2<-inv[3,3] #inv of cov, v2
den1b<-seb1^2+seb3^2
den2b<-2*seb1*seb3
den3b<-pij2/(pii2+pjj2)
denb<-(den1b-(den2b*den3b))^.5
tb<-abs(abs(b1)-abs(b3)) / denb
lambdab<-tb^2
#2 vs 3
pij3<-inv[2,3] #inv of cor between pred of interest 1 vs. 2
pii3<-inv[2,2] #inv of cov, v1
pjj3<-inv[3,3] #inv of cov, v2
den1c<-seb2^2+seb3^2
den2c<-2*seb2*seb3
den3c<-pij3/(pii3+pjj3)
denc<-(den1c-(den2c*den3c))^.5
tc<-abs(abs(b2)-abs(b3)) / denc
lambdac<-tc^2

minusalpha<-1-alpha
Fb<-stats::qf(minusalpha, 1, df)
power12<-round(1-stats::pf(Fb, 1,df,lambdaa),3)
power13<-round(1-stats::pf(Fb, 1,df,lambdab),3)
power23<-round(1-stats::pf(Fb, 1,df,lambdac), 3)
    message("Sample size is ",n)
    message("Power Comparing b1 and b2 = ", power12)
    message("Power Comparing b1 and b3 = ", power13)
    message("Power Comparing b2 and b3 = ", power23)
    result <- data.frame(matrix(ncol = 4))
    colnames(result) <- c( "n","Power b1-b2","Power b1-b3","Power b2-b3")
    result[, 1]<-n
    result[, 2]<-power12
    result[, 3]<-power13
    result[, 4]<-power23
    output<-na.omit(result)
    rownames(output)<- c()
  }
  invisible(output)
  }
