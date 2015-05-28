\name{queryVariant}
\docType{methods}
\alias{queryVariant}
\alias{queryVariant,MyVariant-method}
\alias{queryVariant,missing-method}

\title{
Return  the query result.
}
\description{
This is a wrapper for GET query of "/query?q=<query>" service.
}
\usage{
queryVariant(q, ..., return.as=c("DataFrame", "records", "text"), 
                     myvariant)
}
\arguments{
  \item{q}{
query term(s).
}
  \item{...}{
Commonly queried fields include \code{fields, size} as well as several other fields. 
View available fields by calling \code{?metadata}.
}
  \item{return.as}{
"DataFrame" (default), "records" (list), or "text" (JSON).
}
  \item{myvariant}{
A MyVariant object that describes how to connect to data resources. See \code{\link{MyVariant-class}}. If missing, default object will be used that accesses the main MyVariant.info portal. Default is recommended.
}
}
\value{returns a variant object (DataFrame, list, or JSON text) containing the queried annotations}
\references{
Ref: https://github.com/Network-of-BioThings/myvariant.info/wiki
}

\seealso{
\code{\link{queryVariants}} 
\code{\link{getVariant}} 
\code{\link{getVariants}}
}
\examples{
## return the query result
queryVariant("dbnsfp.variantname:BRCA2")

queryVariant("chr1:1-1000000")
}

