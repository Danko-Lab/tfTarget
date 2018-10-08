\name{mapTF}
\alias{mapTF}
\title{
mapTF
}
\description{
mapTF
}
\usage{
mapTF(tfTar, 
      out.prefix = NULL, 
      distance.cutoff, 
      closest.N,
           gene.pval.cutoff)
}
\arguments{
  \item{tfTar}{tfTarget object obtained by \code\link{diffTXN}}
  \item{out.prefix}{}
  \item{distance.cutoff}{}
  \item{closest.N}{}
}
\value{
  The function adds new items to the input tfTarget object, including:

  \item{TF.TRE.gene.tab }{Data frame indicaiting ...}

}
\examples{
}
\keyword{ tfTarget }
