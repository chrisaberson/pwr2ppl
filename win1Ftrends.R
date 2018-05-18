#'Compute power for a One Factor Within Subjects Trends with up to four levels.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alterative values may be entered by user
#'@param m1 Mean of first time point
#'@param m2 Mean of second time point
#'@param m3 Mean of third time point
#'@param m4 Mean of fourth time point
#'@param s1 Standard deviation of first time point
#'@param s2 Standard deviation of second time point
#'@param s3 Standard deviation of third time point
#'@param s4 Standard deviation of forth time point
#'@param r12 correlation Time 1 and Time 2
#'@param r13 correlation Time 1 and Time 3
#'@param r14 correlation Time 1 and Time 4
#'@param r23 correlation Time 2 and Time 3
#'@param r24 correlation Time 2 and Time 4
#'@param r34 correlation Time 3 and Time 4
#'@param n Sample size for first group
#'@param alpha Type I error (default is .05)
#'@return Power for the One Factor Within Subjects Trends
#'@export

win1Ftrends<-function(m1,m2,m3=NA,m4=NA, s1, s2, s3=NULL,s4=NULL,
                       r12, r13=NULL, r14=NULL, r23=NULL, r24=NULL, r34=NULL,
                       n, alpha=.05)
{
  levels<-NA
  levels[is.na(m4) & is.na(m3)]<-2
  levels[is.na(m4) & !is.na(m3)]<-3
  levels[!is.na(m4)]<-4

  if(levels<3|levels>4){stop("Function requires 3 to 4 levels")}
  if(levels==3){
    var1<-s1^2
    var2<-s2^2
    var3<-s3^2
    cov12<-r12*s1*s2
    cov13<-r13*s1*s3
    cov23<-r23*s2*s3
    out <- mvrnorm(n, mu = c(m1,m2,m3), Sigma = matrix(c(var1,cov12,cov13,
                                                         cov12,var2,cov23,
                                                         cov13, cov23,var3), ncol = 3),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-rename(out, y1 = V1, y2 = V2, y3 = V3)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-gather(out,key="iv",value="dv",-id)
    out$iv<-as.ordered(out$iv)
    options(contrasts=c("contr.helmert", "contr.poly"))

    trends<-aov_car(dv~iv+Error(id/iv),out)
    lin<-testInteractions(trends$lm, custom = list(iv=contr.poly(3,c(1,2,3))[,1]), idata=trends$data[["idata"]])
    lambdaL<-lin$`approx F`
    quad<-testInteractions(trends$lm, custom = list(iv=contr.poly(3,c(1,2,3))[,2]), idata=trends$data[["idata"]])
    lambdaQ<-quad$`approx F`
    minusalpha<-1-alpha
    dfh<-lin$`den Df`*(levels-1)
    dfl<-lin$`den Df`
    FtL<-qf(minusalpha, 1, dfh)
    FtL2<-qf(minusalpha, 1, dfl)
    FtQ<-qf(minusalpha, 1, dfh)
    FtQ2<-qf(minusalpha, 1, dfl)
    powerL.H<-round(1-pf(FtL, 1,dfh,lambdaL),3)
    powerL.L<-round(1-pf(FtL2, 1,dfl,lambdaL),3)
    powerQ.H<-round(1-pf(FtQ, 1,dfh,lambdaQ),3)
    powerQ.L<-round(1-pf(FtQ2, 1,dfl,lambdaQ),3)
    {print(paste("Power Linear Trend for n =",n,"df = ",dfl,"=", powerL.L))}
    {print(paste("Power Linear Trend for n =",n,"df = ",dfh,"=", powerL.H))}
    {print(paste("Power Linear Trend for n =",n,"df = ",dfl,"=", powerQ.L))}
    {print(paste("Power Linear Trend for n =",n,"df = ",dfh,"=", powerQ.H))}}
  if (levels=="4"){
    var1<-s1^2
    var2<-s2^2
    var3<-s3^2
    var4<-s4^2
    cov12<-r12*s1*s2
    cov13<-r13*s1*s3
    cov14<-r14*s1*s4
    cov23<-r23*s2*s3
    cov24<-r24*s2*s4
    cov34<-r34*s3*s4
    out <- mvrnorm(n, mu = c(m1,m2,m3,m4), Sigma = matrix(c(var1,cov12,cov13, cov14,
                                                            cov12,var2,cov23, cov24,
                                                            cov13, cov23,var3, cov34,
                                                            cov14, cov24, cov34, var4), ncol = 4),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-rename(out, y1 = V1, y2 = V2, y3 = V3, y4 = V4)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-gather(out,key="iv",value="dv",-id)
    out$iv<-as.ordered(out$iv)
    options(contrasts=c("contr.helmert", "contr.poly"))
    trends<-aov_car(dv~iv+Error(id/iv),out)
    lin<-testInteractions(trends$lm, custom = list(iv=contr.poly(4,c(1,2,3,4))[,1]), idata=trends$data[["idata"]])
    lambdaL<-lin$`approx F`
    quad<-testInteractions(trends$lm, custom = list(iv=contr.poly(4,c(1,2,3,4))[,2]), idata=trends$data[["idata"]])
    lambdaQ<-quad$`approx F`
    cube<-testInteractions(trends$lm, custom = list(iv=contr.poly(4,c(1,2,3,4))[,3]), idata=trends$data[["idata"]])
    lambdaC<-cube$`approx F`
    minusalpha<-1-alpha
    dfh<-lin$`den Df`*(levels-1)
    dfl<-lin$`den Df`
    FtL<-qf(minusalpha, 1, dfh)
    FtL2<-qf(minusalpha, 1, dfl)
    FtQ<-qf(minusalpha, 1, dfh)
    FtQ2<-qf(minusalpha, 1, dfl)
    FtC<-qf(minusalpha, 1, dfh)
    FtC2<-qf(minusalpha, 1, dfl)
    powerL.H<-round(1-pf(FtL, 1,dfh,lambdaL),3)
    powerL.L<-round(1-pf(FtL2, 1,dfl,lambdaL),3)
    powerQ.H<-round(1-pf(FtQ, 1,dfh,lambdaQ),3)
    powerQ.L<-round(1-pf(FtQ2, 1,dfl,lambdaQ),3)
    powerC.H<-round(1-pf(FtC, 1,dfh,lambdaC),3)
    powerC.L<-round(1-pf(FtC2, 1,dfl,lambdaC),3)
    {print(paste("Power Linear Trend for n =",n,"df = ",dfl,"=", powerL.L))}
    {print(paste("Power Linear Trend for n =",n,"df = ",dfh,"=", powerL.H))}
    {print(paste("Power Quadratic Trend for n =",n,"df = ",dfl,"=", powerQ.L))}
    {print(paste("Power Quadratic Trend for n =",n,"df = ",dfh,"=", powerQ.H))}
    {print(paste("Power Cubic Trend for n =",n,"df = ",dfl,"=", powerC.L))}
    {print(paste("Power Cubic Trend for n =",n,"df = ",dfh,"=", powerC.H))}}
}




