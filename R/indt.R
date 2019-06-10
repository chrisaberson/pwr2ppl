#'Compute power for an Independent Samples t-test
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param m1 Mean of first group
#'@param m2 Mean of second group
#'@param s1 Standard deviation of first group
#'@param s2 Standard deviation of second group
#'@param n1 Sample size for first group
#'@param n2 Sample size for second group
#'@param alpha Type I error (default is .05)
#'@examples
#'indt(m1=22,m2=20,s1=5,s2=5,n1=99,n2=99)
#'indt(m1=1.3, m2=0, s1=4,s2=1,n1=78,n2=234)
#'@return Power for Independent Samples t-test
#'@export
#'
#'
indt<-function(m1=NULL,m2=NULL, s1=NULL,s2=NULL, n1=NULL,n2=NULL, alpha=.05){
  x<-stats::rnorm(n1,m1,s1)
  X<-x
  MEAN<-m1
  SD<-s1
  Z <- (((X - mean(X, na.rm = TRUE))/stats::sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  group<-rep("A1",n1)
  l1<-data.frame(y, group)
  x<-stats::rnorm(n2,m2,s2)
  X<-x
  MEAN<-m2
  SD<-s2
  Z <- (((X - mean(X, na.rm = TRUE))/stats::sd(X, na.rm = TRUE))) * SD
  y<-MEAN + Z
  group<-rep("A2",n2)
  l2<-data.frame(y, group)
  simdat<-rbind(l1,l2)
  anova<-stats::aov(y~group, data=simdat)
  anova<-car::Anova(anova, type="III")
  SSA<-anova[2,1] #column, row
  SSwin<-anova[3,1]
  dfwin<-anova[3,2]
  mswin<-SSwin/dfwin
  dfbg<-anova[2,2]
  eta2<-SSA/(SSA+SSwin)
  f2<-eta2/(1-eta2)
  lambda<-f2*dfwin
  d<-round((m1-m2)/(mswin^.5),3)
  minusalpha<-1-alpha
  Ft<-stats::qf(minusalpha, dfbg, dfwin)
  Power<-1-stats::pf(Ft, dfbg,dfwin,lambda)
  sx1_un <- s1/(n1^.5)
  sx2_un <- s2/(n2^.5)
  sx1x2_un <- ((sx1_un^2)+(sx2_un^2))^.5
  t_un <- (m1-m2)/sx1x2_un
  d_un <- round(t_un * sqrt((n1+n2)/(n1*n2)),3)
  n_harm <- ((2*n1*n2)/(n1+n2))
  var1<-s1^2
  var2<-s2^2
  sat_num <- ((var1/n1)+(var2/n2))*((var1/n1)+(var2/n2))
  sat_denom <- ((((s1^2/n1)^2))/(n1-1)) +((((s2^2/n2)^2))/(n2-1))
  df_un <- round(sat_num/sat_denom,3)
  delta_un <-(d_un*((n_harm/2)^.5))
  lambda_un <-delta_un^2
  Ft_un<-stats::qf(minusalpha, dfbg, df_un)
  Power_unequal<-1-stats::pf(Ft_un, dfbg,df_un,lambda_un)
  pe<-round(Power,3)
  pu<-round(Power_unequal,3)
  message("Equal Variance Power for n1 = ",n1," ,n2 = ",n2, " with d = ", d, " = ", pe)
  message("Unequal Variance Power for n1 = ",n1,"," ," n2 = ",n2, ", with d = ", d_un, " = ", pu)
  result <- data.frame(matrix(ncol = 5))
  colnames(result) <- c("n1", "n2","d","Power Equal Vars","Power Unequal Vars")
  result[, 1]<-n1
  result[, 2]<-n1
  result[, 3]<-d
  result[, 4]<-pe
  result[, 5]<-pu
  output<-na.omit(result)
  rownames(output)<- c()
  invisible(output)
  }
