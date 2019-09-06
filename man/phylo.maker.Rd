\name{phylo.maker}
\alias{phylo.maker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make phylogenetic hypothesis for a list of species based on a backbone phylogeny
}
\description{
This function makes phylogenetic hypotheses for a list of species under three scenarios based on a backbone phylogeny.
}
\usage{
phylo.maker(sp.list, tree = GBOTB.extended, nodes = nodes.info.1, output.sp.list = TRUE, output.tree = FALSE, scenarios = "S3", r = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.list}{An user-provided species list upon which to make the phylogenetic hypotheses.
}
  \item{tree}{The backbone phylogeny. The default is GBOTB.extended.
}
  \item{nodes}{The genus and family nodes information of the backbone phylogeny, the default is nodes.info.1.
}
  \item{output.sp.list}{Whether or not to output the user-provided species list, with a column added to indicate the status of species of pruning, binding or failing to bind to the backbone phylogeny. The default is TRUE.
}
  \item{output.tree}{Whether or not to output the updated backbone phylogeny, within which the user provided species have been incorporated. The default is FALSE.
}
  \item{scenarios}{The scenario under which to make the phylogenetic hypothesis, details can be found in Jin & Qian (2019). The default is S3 (scenario 3).
}
 \item{r}{The number of runs set for S2 (scenario 2), which follows a randomization procedure of inserting tips into the backbone tree. Each run returns a phylogeny built following this procedure. The default is set as 1, can be set to larger integers, such as 100.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of,
  \item{scenario.1 }{The phylogeny of the user specified list of species built under scenario 1.}
  \item{scenario.2 }{A list of phylogenies of the user specified list of species built under scenario 2.}
  \item{scenario.3 }{The phylogeny of the user specified list of species built under scenario 3.}
  \item{species.list }{Optional. The default is TRUE. The user specified list of species with a column at the end indicating of each species status of "prune", "bind" or "fail to bind" to the backbone phylogeny. Specifically, species matched between the species.list and the backbone phylogeny are indicated as "prune"; species absent from the backbone phylogeny, but were successfully binded to the backbone phylogeny are indicated as "bind"; species absent from the backbone phylogeny and failed to be binded to the backbone phylogeny are indicated as  "fail to bind".}
  \item{tree.scenario.1 }{Optional. The backbone phylogeny that has incorporated the species in the user specified species list built under scenario 1. The default is FALSE.}
  \item{tree.scenario.2 }{Optional. A list of the backbone phylogenies that has incorporated the species in the user specified species list built under scenario 2. The default is FALSE.}
  \item{tree.scenario.3 }{Optional. The backbone phylogeny that has incorporated the species in the user specified species list built under scenario 3. The default is FALSE.}
}
\references{
Qian, H. & Jin, Y. (2016) An updated megaphylogeny of plants, a tool for generating plant phylogenies and an analysis of phylogenetic community structure. Journal of Plant Ecology, 9, 233-239.
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
### make the example file
c1 <- c("Carya floridana", "Carya pallida", "Epiprinus siletianus", "Platycarya strobilacea", "Tilia amurensis", "Apodanthes caseariae", "Pilostyles blanchetii")
c2 <- c("Carya", "Carya", "Epiprinus", "Platycarya", "Tilia", "Apodanthes", "Pilostyles")
c3 <- c("Juglandaceae", "Juglandaceae", "Euphorbiaceae", "Juglandaceae", "Malvaceae", "Apodanthaceae", "Apodanthaceae")
example <- data.frame(species = c1, genus = c2, family = c3)

### run the function
result <- phylo.maker(example, scenarios=c("S1","S2","S3"))

### plot the phylogenies with node ages displayed.
library(ape)
par(mfrow = c(1, 3))
plot.phylo(result$scenario.1, cex = 1.5, main = "scenario.1")
nodelabels(round(branching.times(result$scenario.1), 1), cex = 1)
plot.phylo(result$scenario.2[[1]], cex = 1.5, main = "scenario.2")
nodelabels(round(branching.times(result$scenario.2[[1]]), 1), cex = 1)
plot.phylo(result$scenario.3, cex = 1.5, main = "scenario.3")
nodelabels(round(branching.times(result$scenario.3), 1), cex = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
