\name{getAll}
\docType{methods}
\alias{getAll}

\title{
Get all HGVS IDs from Vcf data.frame.
}
\description{
Read in a Vcf data.frame created by \code{\link{getVcf}} to extract all HGVS IDs for querying MyVariant.info. HGVS IDs will be a separate column called "query".
}
\usage{
getAll(vcf.df)
}
\arguments{
  \item{vcf.df}{
Vcf data.frame created by \code{\link{getVcf}}.
}
}
\value{returns a data.frame of your Vcf with all HGVS IDs as new column}

\references{
https://myvariant.info
http://www.hgvs.org/mutnomen/recs-DNA.html
}

\seealso{
\code{\link{getVcf}}
\code{\link{getDels}}  
\code{\link{getSnps}}
\code{\link{getIns}}
\code{\link{getIndels}}
}
\examples{
## return HGVS IDs for all All in a Vcf
all <- getAll(vcf.df)
hgvs <- all$query
}


