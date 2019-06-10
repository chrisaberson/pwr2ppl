#'Compute power for an Chi Square 2x3
#'Takes proportions for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param r1c1 Proportion of overall scores in Row 1, Column 1
#'@param r1c2 Proportion of overall scores in Row 1, Column 2
#'@param r1c3 Proportion of overall scores in Row 1, Column 3
#'@param r2c1 Proportion of overall scores in Row 2, Column 1
#'@param r2c2 Proportion of overall scores in Row 2, Column 2
#'@param r2c3 Proportion of overall scores in Row 2, Column 3
#'@param n Total sample size
#'@param alpha Type I error (default is .05)
#'@examples
#'Chi2X3(r1c1=.25,r1c2=.25,r1c3=.10, r2c1=.10,r2c2=.25,r2c3=.05,n=200)
#'@return Power for 2x3 Chi Square
#'@export
#'
#'

Chi2X3<-function(r1c1, r1c2, r1c3, r2c1, r2c2, r2c3, n, alpha=.05)
{
  df<-2
  po1<-r1c1
  po2<-r1c2
  po3<-r1c3
  po4<-r2c1
  po5<-r2c2
  po6<-r2c3
  pe1<-(r1c1+r1c2+r1c3)*(r1c1+r2c1)
  pe2<-(r1c1+r1c2+r1c3)*(r1c2+r2c2)
  pe3<-(r1c1+r1c2+r1c3)*(r1c3+r2c3)
  pe4<-(r2c1+r2c2+r2c3)*(r1c1+r2c1)
  pe5<-(r2c1+r2c2+r2c3)*(r1c2+r2c2)
  pe6<-(r2c1+r2c2+r2c3)*(r1c3+r2c3)
  lambda<-n*((((po1-pe1)^2)/pe1)+(((po2-pe2)^2)/pe2)+(((po3-pe3)^2)/pe3)+(((po4-pe4)^2)/pe4)+
               (((po5-pe5)^2)/pe5)+(((po6-pe6)^2)/pe6))
  tabled<-stats::qchisq(1-alpha, df=df)
  power<-round(1-stats::pchisq(tabled, df=df, lambda),3)
  sum<-po1+po2+po3+po4+po5+po6
  if(sum!=1.0){stop("Expected proportions must add to 1.0. Check input po values")
  }
  else message("Power for n of ", n, " = ", power)
  result <- data.frame(matrix(ncol = 2))
  colnames(result) <- c( "n","Power")
  result[, 1]<-n
  result[, 2]<-power
  output<-na.omit(result)
  rownames(output)<- c()
  invisible(output)
  }
