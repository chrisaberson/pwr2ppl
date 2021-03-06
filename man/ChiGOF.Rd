% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChiGOF.R
\name{ChiGOF}
\alias{ChiGOF}
\title{Compute power for an Chi Square Goodness of Fit
Takes proportions for up to six group. Alpha is .05 by default, alternative values may be entered by user}
\usage{
ChiGOF(
  groups,
  po1,
  po2,
  po3 = NULL,
  po4 = NULL,
  po5 = NULL,
  po6 = NULL,
  n,
  alpha = 0.05
)
}
\arguments{
\item{groups}{Number of groups}

\item{po1}{Proportion observed Group 1}

\item{po2}{Proportion observed Group 2}

\item{po3}{Proportion observed Group 3}

\item{po4}{Proportion observed Group 4}

\item{po5}{Proportion observed Group 5}

\item{po6}{Proportion observed Group 6}

\item{n}{Total sample size}

\item{alpha}{Type I error (default is .05)}
}
\value{
Power for Chi Square Goodness of Fit
}
\description{
Compute power for an Chi Square Goodness of Fit
Takes proportions for up to six group. Alpha is .05 by default, alternative values may be entered by user
}
\examples{
ChiGOF(po1=.25, po2=.20, po3=.20, po4=.35, groups=4,n=100)
}
