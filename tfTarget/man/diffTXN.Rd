\name{diffTXN}
\alias{diffTXN}
\title{
Select the sigificant TREs from dREG results
}
\description{
Select the siginificant TREs using DEseq package.
}
\usage{
diffTXN(TRE.path, 
      gene.path, 
      bigWig.path, 
      plus.files.query, 
      plus.files.control, 
      minus.files.query, 
      minus.files.control, 
      out.prefix = NULL, 
      ncores = 3)
}
\arguments{
  \item{TRE.path}{File name specified for the TRE results.}
  \item{gene.path}{File name specified gene annotation in bed6 format.}
  \item{bigWig.path}{Folder name indicatingthe the path of the following bigWig file.}
  \item{plus.files.query}{String vector indicating plus bigwig files as CASE.}
  \item{plus.files.control}{String vector indicating plus bigwig files as CONTROL.}
  \item{minus.files.query}{String vector indicating minus bigwig files as CASE.}
  \item{minus.files.control}{String vector indicating minus bigwig files as CONTROL.}
  \item{out.prefix}{String indicating the outpu file prefix.}
  \item{ncores}{Number indicating CPU cores in parallel computing.}
}
\details{
}
\value{
A S3 object with tfTar classname is returned in this function. This object can be passed to follow-uo analysis.
In this step, the following itemms are generated.
  \item{deseq.table.TRE }{BED format including the differentiated TREs}
  \item{deseq.table.gene}{data frame including the genes}
}
\examples{
}
\keyword{ tfTarget }
