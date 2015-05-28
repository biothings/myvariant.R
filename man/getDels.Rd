\name{getDels}
\docType{methods}
\alias{getDels}

\title{
Get deletion HGVS IDs from Vcf data.frame.
}
\description{
Read in a Vcf data.frame created by \code{\link{getVcf}} to extract deletion HGVS IDs for querying MyVariant.info. HGVS IDs will be a separate column called "query".
Only deletions are returned from the original Vcf data.frame.
}
\usage{
getDels(vcf.df)
}
\arguments{
  \item{vcf.df}{
Vcf data.frame created by \code{\link{getVcf}}.
}
}
\value{returns a data.frame of your Vcf with deletion HGVS IDs as new column}

\references{
https://myvariant.info
http://www.hgvs.org/mutnomen/recs-DNA.html
}

\seealso{
\code{\link{getVcf}}
\code{\link{getAll}}  
\code{\link{getSnps}}
\code{\link{getIns}}
\code{\link{getIndels}}
}
\examples{
## return HGVS IDs for all deletions in a Vcf
file.path <- system.file("extdata", "dbsnp_mini.vcf", package="myvariant")
vcf <- getVcf(file.path)
deletions <- getDels(vcf)
hgvs <- deletions$query
}


