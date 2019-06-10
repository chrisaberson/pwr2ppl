#'Compute power for a One Factor ANOVA with three levels.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param m1 Mean of first group
#'@param m2 Mean of second group
#'@param m3 Mean of third group
#'@param s1 Standard deviation of first group
#'@param s2 Standard deviation of second group
#'@param s3 Standard deviation of third group
#'@param n1 Sample size for first group
#'@param n2 Sample size for second group
#'@param n3 Sample size for third group
#'@param alpha Type I error (default is .05)
#'@examples
#'anova1f_3(m1=80, m2=82, m3=82, s1=10, s2=10, s3=10, n1=60, n2=60, n3=60)
#'@return Power for the One Factor ANOVA
#'@export


anova1f_3<-function(m1=NULL,m2=NULL,m3=NULL,s1=NULL,s2=NULL,s3=NULL,n1=NULL,n2=NULL,n3=NULL,alpha=.05){
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
x<-stats::rnorm(n3,m3,s3)
X<-x
MEAN<-m3
SD<-s3
Z <- (((X - mean(X, na.rm = TRUE))/stats::sd(X, na.rm = TRUE))) * SD
y<-MEAN + Z
group<-rep("A3",n3)
l3<-data.frame(y, group)
simdat<-rbind(l1,l2,l3)
anova<-stats::aov(y~group, data=simdat)
anova<-car::Anova(anova, type="III")
SSA<-anova[2,1] #column, row
SSwin<-anova[3,1]
dfwin<-anova[3,2]
dfbg<-anova[2,2]
eta2<-SSA/(SSA+SSwin)
f2<-eta2/(1-eta2)
lambda<-f2*dfwin
minusalpha<-1-alpha
Ft<-stats::qf(minusalpha, dfbg, dfwin)
eta2<-round((eta2),4)
power<-round(1-stats::pf(Ft, dfbg,dfwin,lambda),4)
n<-n1+n2+n3
message("Sample size overall = ",n)
message("Power  = ", power, " for eta-squared = ", eta2)


#Return results silently. Available if written to object
result <- data.frame(matrix(ncol = 6))
colnames(result) <- c("n", "n1","n2", "n3", "Eta-squared", "Power")
result[n, 1]<-n
result[n, 2]<-n1
result[n, 3]<-n2
result[n, 4]<-n3
result[n, 5]<-eta2
result[n, 6]<-power
output<-na.omit(result)
rownames(output)<- c()
invisible(output)

}
