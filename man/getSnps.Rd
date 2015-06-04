\name{getSnps}
\docType{methods}
\alias{getSnps}

\title{
Get Snp HGVS IDs from Vcf data.frame.
}
\description{
Read in a Vcf data.frame created by \code{\link{getVcf}} to extract HGVS IDs for querying MyVariant.info. HGVS IDs will be a separate column called "query".
Only Snps are returned from the original Vcf data.frame.
}
\usage{
getSnps(vcf.df)
}
\arguments{
  \item{vcf.df}{
Vcf data.frame created by \code{\link{getVcf}}.
}
}
\value{returns a data.frame of your Vcf with Snp HGVS IDs as new column}

\references{
https://myvariant.info
http://www.hgvs.org/mutnomen/recs-DNA.html
}

\seealso{
\code{\link{getVcf}}
\code{\link{getAll}}  
\code{\link{getDels}}
\code{\link{getIns}}
}
\examples{
## return HGVS IDs for all Snps in a Vcf
file.path <- system.file("extdata", "dbsnp_mini.vcf", package="myvariant")
vcf <- getVcf(file.path)
snps <- getSnps(vcf)
hgvs <- snps$query
}


