\name{getIns}
\docType{methods}
\alias{getIns}

\title{
Get insertion HGVS IDs from Vcf data.frame.
}
\description{
Read in a Vcf data.frame created by \code{\link{getVcf}} to extract insertion HGVS IDs for querying MyVariant.info. HGVS IDs will be a separate column called "query".
Only insertions are returned from the original Vcf data.frame.
}
\usage{
getIns(vcf.df)
}
\arguments{
  \item{vcf.df}{
Vcf data.frame created by \code{\link{getVcf}}.
}
}
\value{returns a data.frame of your Vcf with insertion HGVS IDs as new column}

\references{
https://myvariant.info
http://www.hgvs.org/mutnomen/recs-DNA.html
}

\seealso{
\code{\link{getVcf}}
\code{\link{getAll}}  
\code{\link{getDels}}
\code{\link{getSnps}}
\code{\link{getIndels}}
}
\examples{
## return HGVS IDs for all insertions in a Vcf
file.path <- system.file("extdata", "dbsnp_mini.vcf", package="myvariant")
vcf <- getVcf(file.path)
insertions <- getIns(vcf)
hgvs <- insertions$query
}


