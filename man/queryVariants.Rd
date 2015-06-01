\name{queryVariants}
\docType{methods}
\alias{queryVariants}
\alias{queryVariants,MyVariant-method}
\alias{queryVariants,missing-method}


\title{
Return the batch query result.
}
\description{
This is a wrapper for POST query of "/query" service.
}
\usage{
queryVariants(qterms, scopes=NULL, ..., 
                      return.as=c("DataFrame", "records", "text"), 
                      myvariant)
}

\arguments{
  \item{qterms}{
A vector or list, or string of comma-separated query terms
}
  \item{scopes}{
Type of types of identifiers, either a list or a comma-separated fields to specify type of
input qterms.
}
  \item{...}{
Commonly queried fields include \code{fields, size} as well as several other fields. 
\code{returnall} returns a list of all related data including duplicated and missing qterms. False by default.
View available fields by calling \code{?metadata}.
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
http://docs.myvariant.info/en/latest/doc/variant_query_service.html#batch-queries-via-post
http://docs.myvariant.info/en/latest/doc/variant_query_service.html#id3
}

\seealso{
\code{\link{queryVariant}}
\code{\link{getVariant}}
\code{\link{getVariants}}
}
\examples{
## return the batch query result
}


