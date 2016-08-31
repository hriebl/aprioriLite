\name{aprioriRules}
\alias{aprioriRules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mining Frequent Associations Rules}
\description{
\code{aprioriRules} mines frequent association rules given a dataset of the
class \code{transactions}, frequent itemsets of the class \code{itemsets}, and
a minimum value for confidence
}
\usage{
aprioriRules(data, allSets, minConf, cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{dataset to work on, should be an object of the class
              \code{transactions}}
  \item{allSets}{itemsets to work on, should be an object of the class
                 \code{itemsets}}
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
sets <- aprioriSets(Groceries, 0.02)
rules <- aprioriRules(Groceries, sets, 0.02)
inspect(rules, "bags")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
