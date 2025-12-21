#'Compute power for a One Factor Within Subjects and One Factor Between ANOVA with up to two by four levels (within).
#'Takes means, sds, and sample sizes for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param m1.1 Mean of first level Within Factor, 1st level Between Factor
#'@param m2.1 Mean of second level Within Factor, 1st level Between Factor
#'@param m3.1 Mean of third level Within Factor, 1st level Between Factor
#'@param m4.1 Mean of fourth level Within Factor, 1st level Between Factor
#'@param m1.2 Mean of first level Within Factor, 2nd level Between Factor
#'@param m2.2 Mean of second level Within Factor, 2nd level Between Factor
#'@param m3.2 Mean of third level Within Factor, 2nd level Between Factor
#'@param m4.2 Mean of fourth level Within Factor, 2nd level Between Factor
#'@param s1.1 Standard deviation of first level Within Factor, 1st level Between Factor
#'@param s2.1 Standard deviation of second level Within Factor, 1st level Between Factor
#'@param s3.1 Standard deviation of third level Within Factor, 1st level Between Factor
#'@param s4.1 Standard deviation of forth level Within Factor, 1st level Between Factor
#'@param s1.2 Standard deviation of first level Within Factor, 2nd level Between Factor
#'@param s2.2 Standard deviation of second level Within Factor, 2nd level Between Factor
#'@param s3.2 Standard deviation of third level Within Factor, 2nd level Between Factor
#'@param s4.2 Standard deviation of forth level Within Factor, 2nd level Between Factor
#'@param r1.2_1 correlation Within Factor Level 1 and Within Factor, Level 2, 1st level Between
#'@param r1.3_1 correlation Within Factor Level 1 and Within Factor, Level 3, 1st level Between
#'@param r1.4_1 correlation Within Factor Level 1 and Within Factor, Level 4, 1st level Between
#'@param r2.3_1 correlation Within Factor Level 1 and Within Factor, Level 3, 1st level Between
#'@param r2.4_1 correlation Within Factor Level 1 and Within Factor, Level 4, 1st level Between
#'@param r3.4_1 correlation Within Factor Level 1 and Within Factor, Level 4, 1st level Between
#'@param r1.2_2 correlation Within Factor Level 1 and Within Factor, Level 2, 2nd level Between
#'@param r1.3_2 correlation Within Factor Level 1 and Within Factor, Level 3, 2nd level Between
#'@param r1.4_2 correlation Within Factor Level 1 and Within Factor, Level 4, 2nd level Between
#'@param r2.3_2 correlation Within Factor Level 1 and Within Factor, Level 3, 2nd level Between
#'@param r2.4_2 correlation Within Factor Level 1 and Within Factor, Level 4, 2nd level Between
#'@param r3.4_2 correlation Within Factor Level 1 and Within Factor, Level 4, 2nd level Between
#'@param r sets same correlations between DVs on all factor levels (seriously, just use this)
#'@param s sets same standard deviation for factor levels (see comment above)
#'@param n for each between group level
#'@param alpha Type I error (default is .05)
#'@examples
#'win1bg1(m1.1 = -.25, m2.1=0, m3.1=0.10, m4.1=.15,
#'m1.2=-.25,m2.2=-.25,m3.2=-.25, m4.2=-.25,
#'s1.1 = .4, s2.1=.5, s3.1=0.6, s4.1=.7,
#'s1.2=.4,s2.2=.5,s3.2=.6, s4.2=.7,n = 50,
#'r1.2_1=.5,r1.3_1=.3,r1.4_1=.15,r2.3_1=.5,r2.4_1=.3,r3.4_1=.5,
#'r1.2_2=.5,r1.3_2=.3,r1.4_2=.15, r2.3_2=.5,r2.4_2=.3,r3.4_2=.5)
#'win1bg1(m1.1 = -.25, m2.1=0, m3.1=0.10, m4.1=.15,
#'m1.2=-.25,m2.2=-.25,m3.2=-.25, m4.2=-.25, s=.4, r = .5, n = 100)
#'@return Power for the One Factor Within Subjects and One Factor Between ANOVA
#'@export
win1bg1<-function(m1.1,m2.1,m3.1=NA,m4.1=NA,m1.2,m2.2,m3.2=NA,m4.2=NA,
                  s1.1=NA,s2.1=NA,s3.1=NA,s4.1=NA,s1.2=NA,s2.2=NA,s3.2=NA,s4.2=NA,
                  r1.2_1=NULL, r1.3_1=NULL, r1.4_1=NULL,
                  r2.3_1=NULL, r2.4_1=NULL,
                  r3.4_1=NULL,
                  r1.2_2=NULL, r1.3_2=NULL, r1.4_2=NULL,
                  r2.3_2=NULL, r2.4_2=NULL,
                  r3.4_2=NULL, r=NULL, s = NULL, n, alpha=.05)

{

  V1<-V2<-V3<-V4<-dv<-ivw<-ivbg<-id<-NULL
  levels<-NA
  levels[is.na(m4.1) & is.na(m3.1)]<-2
  levels[is.na(m4.1) & !is.na(m3.1)]<-3
  levels[!is.na(m4.1)]<-4
  oldoption<-options(contrasts=c("contr.helmert", "contr.poly"))
  oldoption
  on.exit(options(oldoption))

  if (levels=="2"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s1.2<-s;s2.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s1.2^2;var4<-s2.2^2}
    if (!is.null(r)){r1.2_1<-r;r1.2_2<-r}
    cov12<-r1.2_1*s1.1*s2.1
    cov34<-r1.2_2*s1.2*s2.2;
    out1 <- MASS::mvrnorm(n, mu = c(m1.1,m2.1),
                          Sigma = matrix(c(var1,cov12,
                                           cov12,var2), ncol = 2),
                          empirical = TRUE)
    out1<-as.data.frame(out1)
    out1$ivbg<-NA
    out1$ivbg<-1 #identifies group

    out2 <- MASS::mvrnorm(n, mu = c(m1.2,m2.2),
                          Sigma = matrix(c(var3,cov34,
                                           cov34,var4), ncol = 2),
                          empirical = TRUE)
    out2<-as.data.frame(out2)
    out2$ivbg<-NA
    out2$ivbg<-2

    out<-rbind(out1,out2)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, ivbg = ivbg)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="ivw",value="dv",-id, -ivbg)
    out$ivw<-as.factor(out$ivw)
    out$ivbg<-as.factor(out$ivbg)

    model<-afex::aov_ez("id", "dv" , out, within ="ivw", between="ivbg", type=3, detailed=TRUE)
    xx<-summary(model)
    dfA<-xx$univariate.tests[2,2]
    dfB<-xx$univariate.tests[3,2]
    dfAB<-xx$univariate.tests[4,2]
    dfWA<-xx$univariate.tests[2,4]
    dfWB<-xx$univariate.tests[3,4]
    dfWAB<-xx$univariate.tests[4,4]
    SSA<-xx$univariate.tests[2,1]
    SSB<-xx$univariate.tests[3,1]
    SSAB<-xx$univariate.tests[4,1]
    SSWA<-xx$univariate.tests[2,3]
    SSWB<-xx$univariate.tests[3,3]
    SSWAB<-xx$univariate.tests[4,3]
    SSWAB<-model$ANOVA$SSd[4]
    eta2A<-SSA/(SSA+SSWA)
    eta2B<-SSB/(SSB+SSWB)
    eta2AB<-SSAB/(SSAB+SSWAB)
    f2A<-eta2A/(1-eta2A)
    f2B<-eta2B/(1-eta2B)
    f2AB<-eta2AB/(1-eta2AB)
    lambdaA<-f2A*dfWA
    lambdaB<-f2B*dfWB
    lambdaAB<-f2AB*dfWAB
    minusalpha<-1-alpha
    FtA<-stats::qf(minusalpha, dfA, dfWA)
    FtB<-stats::qf(minusalpha, dfB, dfWB)
    FtAB<-stats::qf(minusalpha, dfAB, dfWAB)
    powerA<-round(1-stats::pf(FtA, dfA,dfWA,lambdaA),3)
    powerB<-round(1-stats::pf(FtB, dfB,dfWB,lambdaB),3)
    powerAB<-round(1-stats::pf(FtAB, dfAB,dfWAB,lambdaAB),3)
    eta2A<-round((eta2A),3)
    eta2B<-round((eta2B),3)
    eta2AB<-round((eta2AB),3)
    message("Partial eta-squared Factor A = ", eta2A)
    message("Power Factor A (Between) for n = ",n," is ", powerA)
    message("Partial eta-squared Factor B = ", eta2B)
    message("Power Factor B (Within) for n = ",n," is ", powerB)
    message("Partial eta-squared AxB = ", eta2AB)
    message("Power AxB (Unadjusted) for n = ",n," is ", powerAB)
    message("Both Factors Have 2 levels - There is no adjustment when levels = 2")
    result <- data.frame(matrix(ncol = 7))
    colnames(result) <- c("n", "eta2 A","Power A","eta2 B", "Power B (Within - Unadujsted)",
                          "eta2 AxB","Power AxB")
    result[, 1]<-n
    result[, 2]<-eta2A
    result[, 3]<-powerA
    result[, 4]<-eta2B
    result[, 5]<-powerB
    result[, 6]<-eta2AB
    result[, 7]<-powerAB
    output<-na.omit(result)
    rownames(output)<- c()
  }

  if (levels=="3"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s3.1<-s;s1.2<-s;s2.2<-s;s3.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2;var5<-s^2;var6<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s1.2^2;var5<-s2.2^2; var6<-s3.2^2}
    if (!is.null(r)){r1.2_1<-r;r1.3_1<-r;
    r2.3_1<-r;
    r1.2_2<-r;r1.3_2<-r;
    r2.3_2<-r}
    cov12<-r1.2_1*s1.1*s2.1;cov13<-r1.3_1*s1.1*s3.1;
    cov23<-r2.3_1*s2.1*s3.1;
    cov45<-r1.2_2*s1.2*s2.2;cov46<-r1.3_2*s1.2*s3.2;
    cov56<-r2.3_2*s2.2*s3.2

    out1 <- MASS::mvrnorm(n, mu = c(m1.1,m2.1, m3.1),
                          Sigma = matrix(c(var1,cov12,cov13,
                                           cov12,var2,cov23,
                                           cov13,cov23,var3), ncol = 3),
                          empirical = TRUE)
    out1<-as.data.frame(out1)
    out1$ivbg<-NA
    out1$ivbg<-1 #identifies group

    out2 <- MASS::mvrnorm(n, mu = c(m1.2,m2.2,m3.2),
                          Sigma = matrix(c(var4,cov45,cov46,
                                           cov45,var5,cov56,
                                           cov46,cov56,var6), ncol = 3),
                          empirical = TRUE)
    out2<-as.data.frame(out2)
    out2$ivbg<-NA
    out2$ivbg<-2

    out<-rbind(out1,out2)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2,  y3=V3, ivbg = ivbg)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="ivw",value="dv",-id, -ivbg)
    out$ivw<-as.factor(out$ivw)
    out$ivbg<-as.factor(out$ivbg)
    model<-afex::aov_ez("id", "dv" , out, within ="ivw", between="ivbg", type=3, detailed=TRUE)
    xx<-summary(model)
    dfA<-xx$univariate.tests[2,2]
    dfB<-xx$univariate.tests[3,2]
    dfAB<-xx$univariate.tests[4,2]
    dfWA<-xx$univariate.tests[2,4]
    dfWB<-xx$univariate.tests[3,4]
    dfWAB<-xx$univariate.tests[4,4]
    SSA<-xx$univariate.tests[2,1]
    SSB<-xx$univariate.tests[3,1]
    SSAB<-xx$univariate.tests[4,1]
    SSWA<-xx$univariate.tests[2,3]
    SSWB<-xx$univariate.tests[3,3]
    SSWAB<-xx$univariate.tests[4,3]
    eta2A<-SSA/(SSA+SSWA)
    eta2B<-SSB/(SSB+SSWB)
    eta2AB<-SSAB/(SSAB+SSWAB)
    f2A<-eta2A/(1-eta2A)
    f2B<-eta2B/(1-eta2B)
    f2AB<-eta2AB/(1-eta2AB)
    lambdaA<-f2A*dfWA
    lambdaB<-f2B*dfWB
    lambdaAB<-f2AB*dfWAB
    minusalpha<-1-alpha
    FtA<-stats::qf(minusalpha, dfA, dfWA)
    FtB<-stats::qf(minusalpha, dfB, dfWB)
    FtAB<-stats::qf(minusalpha, dfAB, dfWAB)
    powerA<-round(1-stats::pf(FtA, dfA,dfWA,lambdaA),3)
    powerB<-round(1-stats::pf(FtB, dfB,dfWB,lambdaB),3)
    powerAB<-round(1-stats::pf(FtAB, dfAB,dfWAB,lambdaAB),3)
    ggeB<-xx$pval.adjustments[1,1]
    ggeAB<-xx$pval.adjustments[2,1]
    hfeB<-xx$pval.adjustments[1,3]
    hfeAB<-xx$pval.adjustments[2,3]
    hfeB[hfeB>1]<-1
    hfeAB[hfeAB>1]<-1
    ggdfB<-ggeB*dfB
    ggdfAB<-ggeAB*dfAB
    ggdfWB<-ggeB*dfWB
    ggdfWAB<-ggeAB*dfWAB
    hfdfB<-hfeB*dfB
    hfdfAB<-hfeAB*dfAB
    hfdfWB<-hfeB*dfWB
    hfdfWAB<-hfeAB*dfWAB
    lambdaggB<-f2B*ggdfWB
    lambdaggAB<-f2AB*ggdfWAB
    lambdahfB<-f2B*hfdfWB
    lambdahfAB<-f2AB*hfdfWAB
    FtggB<-stats::qf(minusalpha, ggdfB, ggdfWB)
    FtggAB<-stats::qf(minusalpha, ggdfAB, ggdfWAB)
    FthfB<-stats::qf(minusalpha, hfdfB, hfdfWB)
    FthfAB<-stats::qf(minusalpha, hfdfAB, hfdfWAB)
    powerggB<-round(1-stats::pf(FtggB, ggdfB,ggdfWB,lambdaggB),3)
    powerggAB<-round(1-stats::pf(FtggAB, ggdfAB,ggdfWAB,lambdaggAB),3)
    powerhfB<-round(1-stats::pf(FthfB, hfdfB,hfdfWB,lambdahfB),3)
    powerhfAB<-round(1-stats::pf(FthfAB, hfdfAB,hfdfWAB,lambdaggAB),3)
    eta2A<-round((eta2A),3)
    eta2B<-round((eta2B),3)
    eta2AB<-round((eta2AB),3)
    message("Partial eta-squared Factor A = ", eta2A)
    message("Power Factor A (Between) for n = ",n," is ", powerA)
    message("Partial eta-squared Factor B = ", eta2B)
    message("Power Factor B (Within) for n = ",n," is ", powerB)
    message("Power Factor B H-F Adjusted (Epsilon = ",hfeB, "), for n = " ,n, " is ", powerhfB)
    message("Power Factor B G-G Adjusted (Epsilon = ", ggeB, ") for n = " ,n, " is ", powerggB)
    message("Partial eta-squared Factor AxB = ", eta2AB)
    message("Power AxB (Unadjusted) for n = ",n," is ", powerAB)
    message("Power AxB H-F Adjusted (Epsilon = ",hfeAB ,") for n = ",n, " is ", powerhfAB)
    message("Power AxB G-G Adjusted (Epsilon = ", ggeAB,") for n = ",n, " is ", powerggAB)
    result <- data.frame(matrix(ncol = 15))
    colnames(result) <- c("n", "eta2 A","Power A","eta2 B", "Power B (Within - Unadujsted)", "HF epsilon B",
                          "Power B (HF)","GG Epsilon B","Power B (GG)",
                          "eta2 AxB","Power AxB(Unadjusted)","HF epsilon AxB",
                          "Power AxB(HF)","GG Epsilon AB","Power AxB(GG)")
    result[, 1]<-n
    result[, 2]<-eta2A
    result[, 3]<-powerA
    result[, 4]<-eta2B
    result[, 5]<-powerB
    result[, 6]<-hfeB
    result[, 7]<-powerhfB
    result[, 8]<-ggeB
    result[, 9]<-powerggB
    result[, 10]<-eta2AB
    result[, 11]<-powerAB
    result[, 12]<-hfeAB
    result[, 13]<-powerhfAB
    result[, 14]<-ggeAB
    result[, 15]<-powerggAB
    output<-na.omit(result)
    rownames(output)<- c() }

  if (levels=="4"){
    if (!is.null(s)){
      s1.1<-s; s2.1<-s;s3.1<-s;s4.1<-s;s1.2<-s;s2.2<-s;s3.2<-s;s4.2<-s
      var1<-s^2; var2<-s^2;var3<-s^2;var4<-s^2;var5<-s^2;var6<-s^2;var7<-s^2;var8<-s^2}
    if (is.null(s)){var1<-s1.1^2; var2<-s2.1^2;var3<-s3.1^2;var4<-s4.1^2;var5<-s1.2^2;var6<-s2.2^2;var7<-s3.2^2;var8<-s4.2^2}
    if (!is.null(r)){r1.2_1<-r;r1.3_1<-r;r1.4_1<-r;r2.3_1<-r;r2.4_1<-r;
    r3.4_1<-r;
    r1.2_2<-r;r1.3_2<-r;r1.4_2<-r
    r2.3_2<-r;r2.4_2<-r;
    r3.4_2<-r}
    cov12<-r1.2_1*s1.1*s2.1;cov13<-r1.3_1*s1.1*s3.1;cov14<-r1.4_1*s1.1*s4.1
    cov23<-r2.3_1*s2.1*s3.1;cov24<-r2.4_1*s2.1*s4.1;
    cov34<-r3.4_1*s3.1*s4.1;
    cov56<-r1.2_2*s1.2*s2.2;cov57<-r1.3_2*s1.2*s3.2;cov58<-r1.4_2*s1.2*s4.2
    cov67<-r2.3_2*s2.2*s3.2;cov68<-r2.4_2*s2.2*s4.2
    cov78<-r3.4_2*s3.2*s4.2


    out1 <- MASS::mvrnorm(n, mu = c(m1.1,m2.1, m3.1, m4.1),
                          Sigma = matrix(c(var1,cov12,cov13,cov14,
                                           cov12,var2,cov23,cov24,
                                           cov13,cov23,var3, cov34,
                                           cov14,cov24,cov34,var4), ncol = 4),
                          empirical = TRUE)
    out1<-as.data.frame(out1)
    out1$ivbg<-NA
    out1$ivbg<-1 #identifies group

    out2 <- MASS::mvrnorm(n, mu = c(m1.2,m2.2,m3.2,m4.2),
                          Sigma = matrix(c(var5,cov56,cov57,cov58,
                                           cov56,var6,cov67,cov68,
                                           cov57,cov67,var7,cov78,
                                           cov58,cov68,cov78,var8), ncol = 4),
                          empirical = TRUE)
    out2<-as.data.frame(out2)
    out2$ivbg<-NA
    out2$ivbg<-2

    out<-rbind(out1,out2)
    out<-as.data.frame(out)
    out<-dplyr::rename(out, y1 = V1, y2 = V2, y3=V3, y4=V4,ivbg = ivbg)
    out$id <- rep(1:nrow(out))
    out$id<-as.factor(out$id)
    out<-tidyr::gather(out,key="ivw",value="dv",-id, -ivbg)
    out$ivw<-as.factor(out$ivw)
    out$ivbg<-as.factor(out$ivbg)
    model<-afex::aov_ez("id", "dv" , out, within ="ivw", between="ivbg", type=3, detailed=TRUE)
    xx<-summary(model)
    dfA<-xx$univariate.tests[2,2]
    dfB<-xx$univariate.tests[3,2]
    dfAB<-xx$univariate.tests[4,2]
    dfWA<-xx$univariate.tests[2,4]
    dfWB<-xx$univariate.tests[3,4]
    dfWAB<-xx$univariate.tests[4,4]
    SSA<-xx$univariate.tests[2,1]
    SSB<-xx$univariate.tests[3,1]
    SSAB<-xx$univariate.tests[4,1]
    SSWA<-xx$univariate.tests[2,3]
    SSWB<-xx$univariate.tests[3,3]
    SSWAB<-xx$univariate.tests[4,3]
    eta2A<-SSA/(SSA+SSWA)
    eta2B<-SSB/(SSB+SSWB)
    eta2AB<-SSAB/(SSAB+SSWAB)
    f2A<-eta2A/(1-eta2A)
    f2B<-eta2B/(1-eta2B)
    f2AB<-eta2AB/(1-eta2AB)
    lambdaA<-f2A*dfWA
    lambdaB<-f2B*dfWB
    lambdaAB<-f2AB*dfWAB
    minusalpha<-1-alpha
    FtA<-stats::qf(minusalpha, dfA, dfWA)
    FtB<-stats::qf(minusalpha, dfB, dfWB)
    FtAB<-stats::qf(minusalpha, dfAB, dfWAB)
    powerA<-round(1-stats::pf(FtA, dfA,dfWA,lambdaA),3)
    powerB<-round(1-stats::pf(FtB, dfB,dfWB,lambdaB),3)
    powerAB<-round(1-stats::pf(FtAB, dfAB,dfWAB,lambdaAB),3)
    ggeB<-xx$pval.adjustments[1,1]
    ggeAB<-xx$pval.adjustments[2,1]
    hfeB<-xx$pval.adjustments[1,3]
    hfeAB<-xx$pval.adjustments[2,3]
    hfeB[hfeB>1]<-1
    hfeAB[hfeAB>1]<-1
    ggdfB<-ggeB*dfB
    ggdfAB<-ggeAB*dfAB
    ggdfWB<-ggeB*dfWB
    ggdfWAB<-ggeAB*dfWAB
    hfdfB<-hfeB*dfB
    hfdfAB<-hfeAB*dfAB
    hfdfWB<-hfeB*dfWB
    hfdfWAB<-hfeAB*dfWAB
    lambdaggB<-f2B*ggdfWB
    lambdaggAB<-f2AB*ggdfWAB
    lambdahfB<-f2B*hfdfWB
    lambdahfAB<-f2AB*hfdfWAB
    FtggB<-stats::qf(minusalpha, ggdfB, ggdfWB)
    FtggAB<-stats::qf(minusalpha, ggdfAB, ggdfWAB)
    FthfB<-stats::qf(minusalpha, hfdfB, hfdfWB)
    FthfAB<-stats::qf(minusalpha, hfdfAB, hfdfWAB)
    powerggB<-round(1-stats::pf(FtggB, ggdfB,ggdfWB,lambdaggB),3)
    powerggAB<-round(1-stats::pf(FtggAB, ggdfAB,ggdfWAB,lambdaggAB),3)
    powerhfB<-round(1-stats::pf(FthfB, hfdfB,hfdfWB,lambdahfB),3)
    powerhfAB<-round(1-stats::pf(FthfAB, hfdfAB,hfdfWAB,lambdaggAB),3)
    eta2A<-round((eta2A),3)
    eta2B<-round((eta2B),3)
    eta2AB<-round((eta2AB),3)
    message("Partial eta-squared Factor A = ", eta2A)
    message("Power Factor A (Between) for n = ",n," is ", powerA)
    message("Partial eta-squared Factor B = ", eta2B)
    message("Power Factor B (Within) for n = ",n," is ", powerB)
    message("Power Factor B H-F Adjusted (Epsilon = ",hfeB, "), for n = " ,n, " is ", powerhfB)
    message("Power Factor B G-G Adjusted (Epsilon = ", ggeB, ") for n = " ,n, " is ", powerggB)
    message("Partial eta-squared Factor AxB = ", eta2A)
    message("Power AxB (Unadjusted) for n = ",n," is ", powerAB)
    message("Power AxB H-F Adjusted (Epsilon = ",hfeAB ,") for n = ",n, " is ", powerhfAB)
    message("Power AxB G-G Adjusted (Epsilon = ", ggeAB,") for n = ",n, " is ", powerggAB)
    result <- data.frame(matrix(ncol = 15))
    colnames(result) <- c("n", "eta2 A","Power A","eta2 B", "Power B (Within - Unadujsted)", "HF epsilon B",
                          "Power B (HF)","GG Epsilon B","Power B (GG)",
                          "eta2 AxB","Power AxB(Unadjusted)","HF epsilon AxB",
                          "Power AxB(HF)","GG Epsilon AB","Power AxB(GG)")
    result[, 1]<-n
    result[, 2]<-eta2A
    result[, 3]<-powerA
    result[, 4]<-eta2B
    result[, 5]<-powerB
    result[, 6]<-hfeB
    result[, 7]<-powerhfB
    result[, 8]<-ggeB
    result[, 9]<-powerggB
    result[, 10]<-eta2AB
    result[, 11]<-powerAB
    result[, 12]<-hfeAB
    result[, 13]<-powerhfAB
    result[, 14]<-ggeAB
    result[, 15]<-powerggAB
    output<-na.omit(result)
    rownames(output)<- c()
  }
  invisible(output)}
