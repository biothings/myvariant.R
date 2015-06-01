\name{getVariants}
\docType{methods}
\alias{getVariants}
\alias{getVariants,MyVariant-method}
\alias{getVariants,missing-method}

\title{
Return the list of variant objects for the given list of HGVS ids.
}
\description{
This is a wrapper for POST query of "/variant" service.
}
\usage{
getVariants(hgvsids, fields=NULL, ..., 
          return.as=c("DataFrame", "records", "text"), myvariant)
}

\arguments{
  \item{hgvsids}{
A vector, list, or comm-sep string HGVS ids
}
  \item{fields}{
A vector of fields to return. If fields=="all", all available fields are returned.
}
\item{...}{

}
\item{return.as}{
"DataFrame" (default), "records" (list), "text" (JSON).
}
\item{myvariant}{
A MyVariant object that describes how to connect to data resources. See \code{\link{MyVariant-class}}. If missing, default object will be used that accesses the main MyVariant.info portal. Default is recommended.
}
}

\value{returns a variant object (DataFrame, list, or JSON text) containing the queried annotations}

\references{
http://docs.myvariant.info/en/latest/doc/variant_annotation_service.html#batch-queries-via-post
http://docs.myvariant.info/en/latest/doc/variant_annotation_service.html#id2
}

\seealso{
\code{\link{getVariants}}  
\code{\link{queryVariant}} 
\code{\link{queryVariants}}
}
\examples{
## given a list of HGVS ids
vars <- c('chr1:g.866422C>T',
'chr1:g.876664G>A',
'chr1:g.69635G>C',
'chr1:g.69869T>A',
'chr1:g.881918G>A',
'chr1:g.865625G>A',
'chr1:g.879368C>A',
'chr1:g.889226C>T',
'chr1:g.879492C>G',
'chr1:g.879423T>G',
'chr1:g.881602C>T',
'chr1:g.879115C>G',
'chr1:g.69892T>C',
'chr1:g.879381C>T',
'chr1:g.878330C>G')

## Return the list of variant object for the given list of HGVS ids.
df <- getVariants(vars, fields="dbsnp, wellderly")

}