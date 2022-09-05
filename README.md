prw2ppl Readme
================
Chris Aberson
September 5, 2022

# pwr2ppl

pwr2ppl contains protocols for running a wide range of power analyses.
Analyses range from simple approaches such as t-test and correlations to
multifactor ANOVA (between and within subjects, linear mixed model
approaches) and regression-based approaches (basic multiple regression,
moderated regression, mediation, and logistic regression). Recent additional focus on improving mediation approaches through addition of joint signficance tests and conditional process models. 

These protocols accompany the book [*Applied Power Analysis for the
Behavioral Sciences (2nd
ed.)*](https://www.routledge.com/Applied-Power-Analysis-for-the-Behavioral-Sciences-2nd-Edition-2nd-Edition/Aberson/p/book/9781138044593)

The book gives lots of illustrative examples of using the code but is
not necessary for using the package.

### Prerequisites

I built this under R 4.2.0

## Authors

  - **Chris Aberson** [chrisaberson](https://github.com/chrisaberson)
  
## To install
You will need the devtools package. 

devtools::install_github("chrisaberson/pwr2ppl") 

In some cases, users experience an error message reading as follows:
ERROR: dependencies 'car', 'ez', 'phia', 'afex', 'MBESS', 'lavaan' are not available for package 'pwr2ppl'

Installing car, ez, phia, afex, MBESS, and lavaan manually will rectify this error. 

## License

This project is licensed under the MIT License - see the
[LICENSE.md](LICENSE.md) file for details

[![](https://cranlogs.r-pkg.org/badges/pwr2ppl)](https://cran.r-project.org/package=pwr2ppl)

