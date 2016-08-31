\name{apriori}
\alias{apriori}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mining Frequent Association Rules}
\description{
\code{apriori} mines frequent association rules given a dataset of the class
\code{transactions} and minimum values for both support and confidence
}
\usage{
apriori(data, minSupp, minConf, cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{dataset to work on, should be an object of the class
              \code{transactions}}
  \item{minSupp}{minumum value for support}
  \item{minConf}{minimum value for confidence}
  \item{cl}{optional, \code{snow} cluster for parallel matrix multiplications}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(Groceries)
rules <- apriori(Groceries, 0.02, 0.02)
inspect(rules, "bags")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
