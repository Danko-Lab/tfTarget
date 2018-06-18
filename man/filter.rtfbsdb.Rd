\name{filter.rtfbsdb}
\alias{filter.rtfbsdb}
\title{
filter.rtfbsdb
}
\description{
filter.rtfbsdb}
\usage{
filter.rtfbsdb(tfTar, 
      fdr.cutoff = 0.05, 
      sites.num.cutoff = 10, 
      exp.cutoff = 2, 
      ncores = 1)
}
\arguments{
  \item{tfTar}{tfTarget object obtained by \code\link{searchTFBS} }
  \item{fdr.cutoff}{}
  \item{sites.num.cutoff}{}
  \item{exp.cutoff}{}
  \item{ncores}{Number indicating CPU cores in parallel computing.}
}
\details{
}
\value{
  The function adds new items to the input tfTarget object, including:
  If it is a LIST, use
  \item{TF.TRE.tab }{Description of 'comp1'}
  \item{df.motif.up.all }{Description of 'comp2'}
  \item{df.motif.down.all }{Description of 'comp2'}
}
\examples{

}
\keyword{ tfTarget }