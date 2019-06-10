#'Compute power for an Chi Square 2x2
#'Takes proportions for each group. Alpha is .05 by default, alternative values may be entered by user
#'@param r1c1 Proportion of overall scores in Row 1, Column 1
#'@param r1c2 Proportion of overall scores in Row 1, Column 2
#'@param r2c1 Proportion of overall scores in Row 2, Column 1
#'@param r2c2 Proportion of overall scores in Row 2, Column 2
#'@param n Total sample size
#'@param alpha Type I error (default is .05)
#'@examples
#'Chi2x2(r1c1=.28,r1c2=.22,r2c1=.38,r2c2=.12,n=100)
#'@return Power for 2x2 Chi Square
#'@export
#'
#'

Chi2x2<-function(r1c1, r1c2, r2c1, r2c2, n, alpha=.05)
{
  df<-1 #Defines df
  po1<-r1c1 #Proportion Observed for Row 1, Column 1
  po2<-r1c2 #Proportion Observed for Row 1, Column 2
  po3<-r2c1 #Proportion Observed for Row 2, Column 1
  po4<-r2c2 #Proportion Observed for Row 2, Column 2
  sum<-po1+po2+po3+po4
  pe1<-(r1c1+r1c2)*(r1c1+r2c1) #Proportion Expected for Row 1, Column 1
  pe2<-(r1c1+r1c2)*(r1c2+r2c2)
  pe3<-(r2c1+r2c2)*(r1c1+r2c1)
  pe4<-(r1c1+r1c2)*(r1c2+r2c2)
  lambda<-n*((((po1-pe1)^2)/pe1)+(((po2-pe2)^2)/pe2)+(((po3-pe3)^2)/pe3)+(((po4-pe4)^2)/pe4))
  tabled<-stats::qchisq(1-alpha, df=df)
  power<-round(1-stats::pchisq(tabled, df=df, lambda),3)
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
