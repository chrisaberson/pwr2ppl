#'Compute power for a One Factor Within Subjects Trends with up to four levels.
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
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
#'@examples
#'win1Ftrends(m1=-.25,m2=-.15,m3=-.05,m4=.05,s1=.4,s2=.5,s3=.6,s4=.7,
#'r12=.50, r13=.30, r14=.15, r23=.5, r24=.30, r34=.50, n=25)
#'@import stats
#'@return Power for the One Factor Within Subjects Trends
#'@export

win1Ftrends<-function(m1,m2,m3=NA,m4=NA, s1, s2, s3=NULL,s4=NULL,
                       r12, r13=NULL, r14=NULL, r23=NULL, r24=NULL, r34=NULL,
                       n, alpha=.05)
{
  V1<-V2<-V3<-V4<-dv<-iv<-id<-NULL
  levels<-NA
  levels[is.na(m4) & is.na(m3)]<-2
  levels[is.na(m4) & !is.na(m3)]<-3
  levels[!is.na(m4)]<-4
  options(contrasts=c("contr.helmert", "contr.poly"))
  if(levels<3|levels>4){stop("Function requires 3 to 4 levels")}
  if(levels==3){
    var1<-s1^2
    var2<-s2^2
    var3<-s3^2
    cov12<-r12*s1*s2
    cov13<-r13*s1*s3
    cov23<-r23*s2*s3
    out <- MASS::mvrnorm(n, mu = c(m1,m2,m3), Sigma = matrix(c(var1,cov12,cov13,
                                                         cov12,var2,cov23,
                                                         cov13, cov23,var3), ncol = 3),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3 = V3)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="iv",value="dv",-id)
    out$iv<-as.ordered(out$iv)

    trends<-afex::aov_car(dv~iv+Error(id/iv),out)
    lin<-phia::testInteractions(trends$lm, custom = list(iv=contr.poly(3,c(1,2,3))[,1]), idata=trends$data[["idata"]])
    lambdaL<-lin$`approx F`
    quad<-phia::testInteractions(trends$lm, custom = list(iv=contr.poly(3,c(1,2,3))[,2]), idata=trends$data[["idata"]])
    lambdaQ<-quad$`approx F`
    minusalpha<-1-alpha
    dfh<-lin$`den Df`*(levels-1)
    dfl<-lin$`den Df`
    FtL<-stats::qf(minusalpha, 1, dfh)
    FtL2<-stats::qf(minusalpha, 1, dfl)
    FtQ<-stats::qf(minusalpha, 1, dfh)
    FtQ2<-stats::qf(minusalpha, 1, dfl)
    powerL.H<-round(1-stats::pf(FtL, 1,dfh,lambdaL),3)
    powerL.L<-round(1-stats::pf(FtL2, 1,dfl,lambdaL),3)
    powerQ.H<-round(1-stats::pf(FtQ, 1,dfh,lambdaQ),3)
    powerQ.L<-round(1-stats::pf(FtQ2, 1,dfl,lambdaQ),3)
    message("Power Linear Trend for n = ",n,", df = ",dfl," is ", powerL.L)
    message("Power Linear Trend for n = ",n,", df = ",dfh," is ", powerL.H)
    message("Power Linear Trend for n = ",n,", df = ",dfl," is ", powerQ.L)
    message("Power Linear Trend for n = ",n,", df = ",dfh," is ", powerQ.H)
    result <- data.frame(matrix(ncol = 7))
    colnames(result) <- c("n", "dfl","dfh","Power linear dflow",
                          "Power linear df high","Power quadratic dflow",
                          "Power quadratic df high")
    result[, 1]<-n
    result[, 2]<-dfl
    result[, 3]<-dfh
    result[, 4]<-powerL.L
    result[, 5]<-powerL.H
    result[, 6]<-powerQ.L
    result[, 7]<-powerQ.H
    output<-na.omit(result)
    rownames(output)<- c()
    }
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
    out <- MASS::mvrnorm(n, mu = c(m1,m2,m3,m4), Sigma = matrix(c(var1,cov12,cov13, cov14,
                                                            cov12,var2,cov23, cov24,
                                                            cov13, cov23,var3, cov34,
                                                            cov14, cov24, cov34, var4), ncol = 4),
                   empirical = TRUE)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3 = V3, y4 = V4)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="iv",value="dv",-id)
    out$iv<-as.ordered(out$iv)
    trends<-afex::aov_car(dv~iv+Error(id/iv),out)
    lin<-phia::testInteractions(trends$lm, custom = list(iv=contr.poly(4,c(1,2,3,4))[,1]), idata=trends$data[["idata"]])
    lambdaL<-lin$`approx F`
    quad<-phia::testInteractions(trends$lm, custom = list(iv=contr.poly(4,c(1,2,3,4))[,2]), idata=trends$data[["idata"]])
    lambdaQ<-quad$`approx F`
    cube<-phia::testInteractions(trends$lm, custom = list(iv=contr.poly(4,c(1,2,3,4))[,3]), idata=trends$data[["idata"]])
    lambdaC<-cube$`approx F`
    minusalpha<-1-alpha
    dfh<-lin$`den Df`*(levels-1)
    dfl<-lin$`den Df`
    FtL<-stats::qf(minusalpha, 1, dfh)
    FtL2<-stats::qf(minusalpha, 1, dfl)
    FtQ<-stats::qf(minusalpha, 1, dfh)
    FtQ2<-stats::qf(minusalpha, 1, dfl)
    FtC<-stats::qf(minusalpha, 1, dfh)
    FtC2<-stats::qf(minusalpha, 1, dfl)
    powerL.H<-round(1-stats::pf(FtL, 1,dfh,lambdaL),3)
    powerL.L<-round(1-stats::pf(FtL2, 1,dfl,lambdaL),3)
    powerQ.H<-round(1-stats::pf(FtQ, 1,dfh,lambdaQ),3)
    powerQ.L<-round(1-stats::pf(FtQ2, 1,dfl,lambdaQ),3)
    powerC.H<-round(1-stats::pf(FtC, 1,dfh,lambdaC),3)
    powerC.L<-round(1-stats::pf(FtC2, 1,dfl,lambdaC),3)
    message("Power Linear Trend for n = ",n,", df = ",dfl," is ", powerL.L)
    message("Power Linear Trend for n = ",n,", df = ",dfh," is ", powerL.H)
    message("Power Quadratic Trend for n = ",n,", df = ",dfl," is ", powerQ.L)
    message("Power Quadratic Trend for n = ",n,", df = ",dfh," is ", powerQ.H)
    message("Power Cubic Trend for n = ",n,", df = ",dfl," is ", powerC.L)
    message("Power Cubic Trend for n = ",n,", df = ",dfh," is ", powerC.H)
    result <- data.frame(matrix(ncol = 9))
    colnames(result) <- c("n", "dfl","dfh","Power linear dflow",
                          "Power linear df high","Power quadratic dflow",
                          "Power quadratic df high","Power cubic dflow",
                          "Power cubic df high")
    result[, 1]<-n
    result[, 2]<-dfl
    result[, 3]<-dfh
    result[, 4]<-powerL.L
    result[, 5]<-powerL.H
    result[, 6]<-powerQ.L
    result[, 7]<-powerQ.H
    result[, 8]<-powerC.L
    result[, 9]<-powerC.H
    output<-na.omit(result)
    rownames(output)<- c()
  }
  invisible(output)
    }




