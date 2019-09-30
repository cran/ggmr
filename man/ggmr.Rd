\name{ggmr}
\alias{ggmr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solves the Generalized Gauss Markov Regression model
}
\description{
Fits the linear model using covariance matrices on the predictor, the response and covariance matrix between predictor and response, according to ISO/TS 28037 (2010).
}
\usage{
ggmr(x, y, Ux = diag(0, length(x)), 
			Uy = diag(1, length(x)), 
			Uxy = diag(0, length(x)),
			subset = rep(TRUE, length(x)), 
			tol = sqrt(.Machine$double.eps), max.iter = 100, alpha = 0.05,
			coef.H0 = c(0, 1))}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric vector, the predictor values}
  \item{y}{numeric vector, the response values}
  \item{Ux}{numeric matrix, the variance matrix of the predictor}
  \item{Uy}{numeric matrix, the variance matrix of the response}
  \item{Uxy}{numeric matrix, the covariance matrix between predictor and the response}
  \item{subset}{a logical vector or a numeric vector with the position to be considered}
  \item{tol}{numeric, the maximum allowed error tolerance, tolerance is relative}
  \item{max.iter}{integer, the maximum number of allowed iterations}
  \item{alpha}{numeric, the significance level used on testing H0}
  \item{coef.H0}{the coeffients for hypothesis testing purposes}
}
\value{
	a list with the following elements
  \item{coefficients}{estimated coefficients} 
  \item{cov}{covariance matrix of the estimated coefficients}
  \item{xi}{estimated latent unobservable variables}
  \item{chisq.validation}{chi-squared statistic for model validation}
  \item{chisq.ht}{chi-squared statistic of the observed values for the hypothesis testing}
  \item{chisq.cri}{chi-squared critical value}
  \item{p.value}{probability of observing a validation statistic equal or larger then the sampled just by chance}
  \item{curr.iter}{current number of iterations used}
  \item{curr.tol}{current relative tolerance} 
}
\references{
ISO/TS 28037 (2010). \emph{Determination and Use of straight-line calibration functions} \url{https://www.iso.org/standard/44473.html}
}
\author{
Hugo Gasca-Aragon

Maintainer: Hugo Gasca-Aragon <hugo_gasca_aragon@hotmail.com>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{lm}}, \code{\link[dwlm]{dwlm}}
}
\examples{

require(MASS)

# Example ISO 28037 (2010) Section 6. table 6
d<- data.frame(
 x=c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
 y=c(3.2, 4.3, 7.6, 8.6, 11.7, 12.8),
 uy=c(0.5, 0.5, 0.5, 1.0, 1.0, 1.0)
)
# estimates
ggmr.res <- ggmr(d$x, d$y, Uy=diag(d$uy^2), coef.H0=c(0, 2), tol = 1e-10)
ggmr.res$coefficients
sqrt(diag(ggmr.res$cov))
ggmr.res$cov[1, 2]
ggmr.res$chisq.validation
ggmr.res$chisq.cri
# reference values
# coefficients = c(0.885, 2.057)
# se = c(0.530, 0.178)
# cov = -0.082
# validation.stat = 4.131
# critical.value = 9.488

# lm() estimates the coefficients correctly but 
# fails to reproduce the standard errors
summary(lm(y~x, data=d, weights=1/d$uy^2)) 
# coefficients = c(0.8852, 2.0570)
# se = c(0.5383, 0.1808)

# Example ISO 28037 (2010) Section 7. table 10
d <- data.frame(
 x = c(1.2, 1.9, 2.9, 4.0, 4.7, 5.9),
 y = c(3.4, 4.4, 7.2, 8.5, 10.8, 13.5)
)
Ux = diag(c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2))^2
Uy = diag(c(0.2, 0.2, 0.2, 0.4, 0.4, 0.4))^2

# estimates
ggmr.res <- ggmr(d$x, d$y, Ux, Uy, coef.H0=c(0, 2), tol = 1e-10)
ggmr.res$coefficients
sqrt(diag(ggmr.res$cov))
ggmr.res$cov[1, 2]
ggmr.res$chisq.validation
ggmr.res$chisq.cri
# reference values
# coefficients = c(0.5788, 2.1597)
# se = c(0.4764, 0.1355)
# cov = -0.0577
# validation.stat = 2.743
# critical.value = 9.488

# Example ISO 28037 (2010) Section 10. table 25
d<- data.frame(
 x=c(50.4, 99.0, 149.9, 200.4, 248.5, 299.7, 349.1),
 y=c(52.3, 97.8, 149.7, 200.1, 250.4, 300.9, 349.2)
)

Ux<- matrix(c(
 0.50, 0.00, 0.25, 0.00, 0.25, 0.00, 0.25,
 0.00, 1.25, 1.00, 0.00, 0.00, 1.00, 1.00,
 0.25, 1.00, 1.50, 0.00, 0.25, 1.00, 1.25,
 0.00, 0.00, 0.00, 1.25, 1.00, 1.00, 1.00,
 0.25, 0.00, 0.25, 1.00, 1.50, 1.00, 1.25,
 0.00, 1.00, 1.00, 1.00, 1.00, 2.25, 2.00,
 0.25, 1.00, 1.25, 1.00, 1.25, 2.00, 2.50
), 7, 7)

Uy<- matrix(1.00, 7, 7) + diag(4.00, 7)

Uxy<- matrix(0, 7, 7)

# estimates
ggmr.res<- ggmr(d$x, d$y, Ux, Uy, Uxy)
ggmr.res$coefficients
sqrt(diag(ggmr.res$cov))
ggmr.res$cov[1, 2]
ggmr.res$chisq.validation
ggmr.res$chisq.cri
# reference values
# coefficients = c(0.3424, 1.0012)
# se = c(2.0569, 0.0090)
# cov = -0.0129
# validation.stat = 1.772
# critical.value = 11.070

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ggmr }
% __ONLY ONE__ keyword per line
