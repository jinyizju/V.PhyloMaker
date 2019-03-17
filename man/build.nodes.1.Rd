\name{build.nodes.1}
\alias{build.nodes.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the genus- and family-level nodes information in a phylogeny
}
\description{
This function extracts the genus- and family-level nodes information in a phylogeny.
}
\usage{
build.nodes.1(tree, tips)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
The phylogeny.
}
  \item{tips}{
The tips information, includes all the species (i.e., tip labels) in the phylogeny, as well as their genus and family assignments. This file should include three columns named as "species", "genus" and "family".
}
}
\details{
Extracts the genus- and family-level basal and root node labels, node ages, etc., in a phylogeny, to be used in functions such as 'phylo.maker' and 'bind.relative'. Specifically, the basal node of a clade is defined as the most recent common ancestor(MRCA) of two most distantly related tips of that clade in the phylogeny, and the root node is the node rightly above the basal node. In the extreme case when one clade owns only one tip in the phylogeny, the basal node is the same as the root node, which is the node the tip driectly roots from. If a clade is distributed in more than one clusters(i.e. polyphyletic), this function extracts the nodes information of the largest cluster to represent the clade.
}
\value{
The data frame of nodes information, including 10 columns.
  \item{level }{The level of the clade, be either F or G. F = family, G = genus.}
  \item{family }{The family name.}
  \item{genus }{The genus name.}
  \item{rn }{The root node label.}
  \item{rn.bl }{The root node age.}
  \item{bn }{The basal node label.}
  \item{bn.bl }{The basal node age.}
  \item{gen.n }{The number of genera that the clade includes.}
  \item{sp.n }{The number of species that the clade includes.}
  \item{taxa }{The species name, in case the clade includes only one species.}
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
# make the tree and its tips information
phylogeny <- rcoal(10)  # generate a phylogeny with 10 tips
phylogeny$node.label <- paste("N",1:phylogeny$Nnode,sep="")  # label the nodes
phylogeny$tip.label <- paste(rep(rep(c("Genus1","Genus2"),2),4:1), phylogeny$tip.label, sep="_")   # rename the tips
tips <- data.frame(species = phylogeny$tip.label,genus=rep(rep(c("Genus1","Genus2"),2),4:1),family=rep(rep(c("Family1","Family2"),2),4:1)) # generate the tips information

# plot the tree
plot(phylogeny, show.node.label=TRUE, tip.color=rep(rep(c("black","red"),2),4:1),edge.width=1) # plot the phylogeny

# generate the nodes information of the tree
nodes.inf<-build.nodes.1(phylogeny, tips)   # build the nodes information file for the phylogeny
nodes.inf   # show the nodes.infomation file
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
