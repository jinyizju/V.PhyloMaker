\name{bind.relative}
\alias{bind.relative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bind tips to their designated relatives
}
\description{
This function binds a list of species to their designated species- and genus-level relatives in the phylogeny.
}
\usage{
bind.relative(sp.list, tree = GBOTB.extended, nodes = nodes.info.1, output.sp.list = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.list}{
The user-provided species list.
}
  \item{tree}{
The backbone phylogeny. The default is GBOTB.extended.
}
  \item{nodes}{
The nodes information of the backbone phylogeny. The default is nodes.info.1
}
  \item{output.sp.list}{
Whether or not to output the species list, with each species' status of binded to its designated relative or not. The default is TRUE.
}
}
\details{
}
\value{

  \item{phylo }{The updated phylogeny, with species binded to the designated relative}
  \item{species.list }{The species list with the information of each species status of binded to its designated relative or not.}
  \item{nodes.info }{The nodes information of the updated phylogeny.}

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
### make the example file
c1 <- c("Carya floridana", "Carya pallida", "Epiprinus siletianus", "Platycarya strobilacea", "Tilia amurensis", "Apodanthes caseariae", "Pilostyles blanchetii")
c2 <- c("Carya", "Carya", "Epiprinus", "Platycarya", "Tilia", "Apodanthes", "Pilostyles")
c3 <- c("Juglandaceae", "Juglandaceae", "Euphorbiaceae", "Juglandaceae", "Malvaceae", "Apodanthaceae", "Apodanthaceae")
c4 <- c("", "Carya floridana", "", "", "", "", "")
c5 <- c("", "", "", "Cyclocarya", "", "", "")
example <- data.frame(species = c1, genus = c2, family = c3, species.relative = c4, genus.relative = c5)

### run the function with phylo.maker.
result <- bind.relative(example)
result1 <- phylo.maker(sp.list = result$species.list, tree = result$phylo, nodes = result$nodes.info, scenarios = "S3")

### run phylo.maker only.
result2 <- phylo.maker(sp.list = example, scenarios="S3")

### compare the phylogeny with and without incorporting bind.relative. node age is displayed.
par(mfrow=c(1,2))
plot.phylo(result1$scenario.3, cex =1.5, main = "bind.relative + phylo.maker")
nodelabels(round(branching.times(result1$scenario.3), 1), cex = 1)
plot.phylo(result2$scenario.3, cex = 1.5, main = "phylo.maker")
nodelabels(round(branching.times(result2$scenario.3), 1), cex = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
