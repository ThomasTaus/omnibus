\name{hc.test}
\alias{hc.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Higher Criticism
}
\description{
\code{hc-test} computes p-values under the global null hypothesis using Higher Criticism.
}
\usage{
hc.test(p, tuning = c("half", "halfmin", "all"), N.sim = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
numeric vector or matrix of input p-values. In the case of a matrix, rows correspond to individual hypotheses.
}
  \item{tuning}{
character string specifying the tuning parameter (see 'Details'). Possible values are \code{'half'}, \code{'halfmin'} or \code{'all'}.
}
  \item{N.sim}{
numeric indicating the number of simulations that should be performed under H0.
}
}
\details{
After sorting p-values in ascending order the Higher Criticism statistic is computed depending on the tuning parameter defined by \code{tuning}:
\describe{
  \item{half}{only the most extreme half of p-values are used.}
  \item{halfmin}{similar to \code{'half'} but also excluding the most extreme p-value.}
  \item{all}{all p-values are used.}
}
}
\value{
A numeric vector of p-values under the global null hypothesis, with length of 1 if \code{p} is a vector or length equal \code{ncol(p)} if \code{p} is a matrix.
}
\references{
Tukey J.W.: T13: N the higher criticism. Course nodes, \emph{Statistics 411, Princeton University} \bold{1976}.

Donoho D. and Jin J.: Higher criticism for detecting sparse heterogeneous mixtures, \emph{Ann. Stat.} \bold{2004}, 32:962-994.

Donoh D. and Jin J.: Higher criticism for large-scale inference: especially for rare and weak effects, \emph{Stat. Sci.} \bold{2015}, 30:1-25.
}
\author{
Thomas Taus and Andreas Futschik
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{omnibus.test}}.
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
