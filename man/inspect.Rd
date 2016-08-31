\name{inspect}
\alias{inspect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Inspecting Frequent Itemsets or Association Rules}
\description{
\code{inspect} allows the user to find out more about frequent itemsets or
association rules, usually returned from \code{apriori}, \code{aprioriSets},
or \code{aprioriRules}
}
\usage{
inspect(data, search, minSupp, minConf, minLift, minConv, sortBy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{frequent itemsets or association rules to work on, should be an
              object of the class \code{itemsets} or \code{rules}}
  \item{search}{search string for item labels}
  \item{minSupp}{minimum value for support}
  \item{minConf}{minimum value for confidence, not applicable for itemsets}
  \item{minLift}{minimum value for lift, not applicable for itemsets}
  \item{minConv}{minimum value for conviction, not applicable for itemsets}
  \item{sortBy}{sort itemsets or rules by quality measure, can be
                \code{"support"} for itemsets and \code{"support"},
                \code{"confidence"}, \code{"lift"}, or \code{"conviction"}
                for rules}
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
inspect(rules, "bags", minSupp = 0.024)
inspect(rules, "bags", minSupp = 0.024, sortBy = "lift")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
