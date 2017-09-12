\name{omnibus.test}
\alias{omnibus.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Omnibus Test for the Global Null Hypothesis
}
\description{
\code{omnibus.test} tests the global null hypothesis. Rejection can be interpreted as evidence for at least one alternative hypothesis to be true.
}
\usage{
omnibus.test(p, method = c("z", "p", "log.p"), N.sim = 10000, approximate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
numeric vector or matrix of input p-values. In the case of a matrix, rows correspond to individual hypotheses.
}
  \item{method}{
character string naming a pre-defined transformation of input p-values (see 'Details'). Alternatively, a custom function can be specified that takes a numeric vector as input and returns a vector of transformed values.
}
  \item{N.sim}{
numeric indicating the number of simulations that should be performed under H0.
}
  \item{approximate}{
logical specifying whether a beta-approximation should be applied to allow for extreme p-values, or not (default).
}
}
\details{
Before summing up p-values they are transformed according to \code{method}. The following pre-defined options are available:
\describe{
  \item{z}{p-values are transformed to z-scores.}
  \item{p}{\code{1-p}.}
  \item{log.p}{\code{-log(p)}.}
}
}
\value{
A numeric vector of p-values under the global null hypothesis, with length of 1 if \code{p} is a vector or length equal \code{ncol(p)} if \code{p} is a matrix.
}
\references{
Futschik A., Taus T. and Zehetmayer S.: An omnibus test for the global null hypothesis, \emph{arXiv} \bold{2017}, 1709.00960.
}
\author{
Thomas Taus and Andreas Futschik
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hc.test}}.
}
\examples{
# perform Omnibus-test for 5 p-values under H0
p <- runif(5)
omnibus.test(p=p)

# perform 10,000 Omnibus-tests for p-values under H0
n <- 10000
p <- matrix(runif(n*5), nrow=5)
p.global <- omnibus.test(p=p)
hist(p.global)

# test 3 hypotheses with and without effekt
p1 <- t.test(rnorm(15, mean=1), rnorm(15, mean=1))$p.value # mean identical
p2 <- t.test(rnorm(15, mean=2), rnorm(15, mean=2))$p.value # mean identical
p3 <- t.test(rnorm(15, mean=1), rnorm(15, mean=2))$p.value # mean different
# perform Omnibus-test of the null hypothesis that none of the three tests is significant
omnibus.test(p=c(p1, p2, p3))

# repeat Omnibus-test using another predefined method to transform input p-values
omnibus.test(p=c(p1, p2, p3), method="log.p")

# repeat Omnibus-test using a custom transformation of input p-values
transformation <- function(x) {
  x^-2
}
omnibus.test(p=c(p1, p2, p3), method=transformation)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
