\name{nodes.info.1}
\alias{nodes.info.1}
\docType{data}
\title{
The nodes information of 'GBOTB.extended' extracted by 'build.nodes.1'
}
\description{
The genus- and family-level nodes information of 'GBOTB.extended' extracted by the function 'build.nodes.1'.
}
\usage{}
\format{
  A data frame with 11083 rows on the following 10 variables.
  \describe{
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
}
\details{
This data frame is the output of 'build.nodes.1', with the input of 'GBOTB.extended' and 'tips.info'.


}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
}
\keyword{datasets}
