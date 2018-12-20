#'Compute power for a One Factor Between Subjects ANOVA with four levels
#'Takes means, sds, and sample sizes for each group
#'@param m1 Mean of first group
#'@param m2 Mean of second group
#'@param m3 Mean of third group
#'@param m4 Mean of fourth group
#'@param s1 Standard deviation of first group
#'@param s2 Standard deviation of second group
#'@param s3 Standard deviation of third group
#'@param s4 Standard deviation of forth group
#'@param n1 Sample size for first group
#'@param n2 Sample size for second group
#'@param n3 Sample size for third group
#'@param n4 Sample size for fourth grou
#'@param alpha Type I error (default is .05)
#'@return Power for the One Factor Between Subjects ANOVA
#'@export

anova1f_4<-function(m1=NULL,m2=NULL,m3=NULL,m4=NULL, s1=NULL,s2=NULL,s3=NULL,s4=NULL,
                         n1=NULL,n2=NULL,n3=NULL,n4=NULL, alpha=.05){
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
x<-stats::rnorm(n4,m4,s4)
X<-x
MEAN<-m4
SD<-s4
Z <- (((X - mean(X, na.rm = TRUE))/stats::sd(X, na.rm = TRUE))) * SD
y<-MEAN + Z
group<-rep("A4",n4)
l4<-data.frame(y, group)
simdat<-rbind(l1,l2,l3,l4)
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
power<-round(1-stats::pf(Ft, dfbg,dfwin,lambda),3)
#list(Power = power)
{print(paste("Power =", power))}
}
