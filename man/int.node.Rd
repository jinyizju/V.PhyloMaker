\name{int.node}
\alias{int.node}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bind a tip above a node
}
\description{
This functions binds a tip to an interior branch above a node, and creates a node.
}
\usage{
int.node(phylogeny, location.node, tip.label, node.label = NULL, position = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phylogeny}{
The phylogeny.
}
  \item{location.node}{
The number or label of the node above which the tip is to bind.
}
  \item{tip.label}{
The label of the tip to bind.
}
  \item{node.label}{
The label of the node created when binding the tip. The default is NULL.
}
  \item{position}{
The position on the interior branch to bind the new tip, the value must be > 0 and < 1, meaning the proportion of the interior branch length starts from location.node. The default is 0.5.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{phylogeny }{The updated phylogeny, with the tip binded.}
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
phylo0 <- rcoal(10)    # generate a phylogeny with 10 tips
phylo <- int.node(phylogeny = phylo0, location.node = sample(phylo0$edge[, 1], 1), tip.label = "New Tip", node.label = NULL, position = 0.5)    # bind a "New Tip" to the phylogeny to a branch above a random node.
par(mfrow = c(1, 2))
plot(phylo0)    # plot the phylogeny.
plot(phylo)    # plot the phylogeny with the "New Tip" binded.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
