\name{getIndels}
\docType{methods}
\alias{getIndels}

\title{
Get deletion/insertion HGVS IDs from Vcf data.frame.
}
\description{

Read in a Vcf data.frame created by \code{\link{getVcf}} to extract deletion/insertion HGVS IDs for querying MyVariant.info. HGVS IDs will be a separate column called "query". Only deletion/insertions are returned from the original Vcf data.frame.
}
\usage{
getIndels(vcf.df)
}
\arguments{
  \item{vcf.df}{
Vcf data.frame created by \code{\link{getVcf}}.
}
}
\value{returns a data.frame of your Vcf with deletion/insertion (delins) HGVS IDs as new column}

\references{
Ref: https://myvariant.info
http://www.hgvs.org/mutnomen/recs-DNA.html
}

\seealso{
\code{\link{getVcf}}
\code{\link{getAll}}  
\code{\link{getDels}}
\code{\link{getIns}}
\code{\link{getSnps}}
}
\examples{
## return HGVS IDs for all Indels in a Vcf
indels <- getIndels(vcf.df)
hgvs <- indels$query
}


