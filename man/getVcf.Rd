\name{getVcf}
\docType{methods}
\alias{getVcf}

\title{
Read in a Vcf file as a data.frame object.
}
\description{
This is a simple function to read in your Vcf file only for the purpose of passing it to readAll, readSnps, readDels, readIns, or readIndels.
The five functions will extract information to creat proper HGVS IDs for querying MyVariant.info.
}
\usage{
getVcf(file.path)
}
\arguments{
  \item{file.path}{
path to Vcf file.
}
}
\value{returns a data.frame of your Vcf}

\references{
Ref: https://myvariant.info
}

\seealso{
\code{\link{getAll}} 
\code{\link{getSnps}} 
\code{\link{getDels}}
\code{\link{getIns}}
}
\examples{
## return Vcf data.frame
file.path <- system.file("extdata", "dbsnp_mini.vcf", package="myvariant")
vcf.df <- getVcf(file.path)
}


