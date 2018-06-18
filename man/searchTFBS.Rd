\name{searchTFBS}
\alias{searchTFBS}
\title{
Search TFBS 
}
\description{
Search TFBS 
}
\usage{
searchTFBS(tfTar, 
      tfs, 
      file.twoBit, 
      pval.cutoff.up = 0.01, 
      pval.cutoff.down = 0.1, 
      half.size = 150, 
      mTH = 7, 
      min.size = 150, 
      run.repeats = 2, 
      ncores = 1)
}
\arguments{
  \item{tfTar}{tfTarget object obtained by \code\link{diffTXN}}
  \item{tfs}{RDATA file including the TFBS identified by RtfbsDB package.}
  \item{file.twoBit}{File name indicating the binary data of sequence.}
  \item{pval.cutoff.up}{Float indicating.}
  \item{pval.cutoff.down}{Float indicating.}
  \item{half.size}{Number indicating.}
  \item{mTH}{Float indicating the threshold used in RtfbsDB methods.}
  \item{min.size}{Number indicating.}
  \item{run.repeats}{Number indicating.}
  \item{ncores}{Number indicating CPU cores in parallel computing.}
}
\details{
}
\value{
  The function adds new items to the input tfTarget object, including:
  \item{enh.up.bed }{Description of 'comp1'}
  \item{enh.down.bed }{Description of 'comp2'}
  \item{enh.unc.bed }{Description of 'comp2'}
  \item{motif.list.up }{Description of 'comp2'}
  \item{motif.list.down }{Description of 'comp2'}
}
\examples{
}
\keyword{ tfTarget }
