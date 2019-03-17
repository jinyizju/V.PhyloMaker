\name{at.node}
\alias{at.node}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bind a tip to a node
}
\description{
This function binds a tip to a node and creates a polytomy
}
\usage{
at.node(phylogeny, location.node, tip.label)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{
The phylogeny.
}
  \item{location.node}{
The number or label of the node to bind.
}
  \item{tip.label}{
The label of the tip.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 Returns the phylogeny with the tip binded to the node.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yi Jin
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(ape)
phylo0 <- rcoal(10)  # generate a phylogeny with 10 tips
phylo <- at.node(phylogeny = phylo0, location.node = sample(phylo0$edge[, 1], 1), tip.label = "New Tip")  # bind a "New Tip" to the phylogeny at a random node.
par(mfrow = c(1, 2))
plot(phylo0)        # plot the phylogeny.
plot(phylo)          # plot the phylogeny with the "New Tip" binded.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
