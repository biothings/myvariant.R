\name{getVariant}
\docType{methods}
\alias{getVariant}
\alias{getVariant,MyVariant-method}
\alias{getVariant,missing-method}

\title{
Return the variant object for the given HGVS id.
}
\description{
This is a wrapper for GET query of "/variant/<hgvsid>" service.
}
\usage{
getVariant(hgvsid, fields=NULL, 
        ..., return.as=c("records", "text"), myvariant)
}
\arguments{
  \item{hgvsid}{
HGVS id
}
  \item{fields}{
Fields to return, a list of a comma-sep string. If fields=="all", all available fields are returned.
}
\item{...}{

}
\item{return.as}{
"records"" (list), "text" (JSON).
}
\item{myvariant}{
A MyVariant object that describes how to connect to data resources. See \code{\link{MyVariant-class}}. If missing, default object will be used that accesses the main MyVariant.info portal. Default is recommended.
}
}
\value{returns a variant object containing the queried annotations}

\references{
Ref: https://github.com/Network-of-BioThings/myvariant.info/wiki
}

\seealso{
\code{\link{getVariants}} 
\code{\link{queryVariant}} 
\code{\link{queryVariants}}
}
\examples{
## return the variant object for the given HGVS id
getVariant("chr7:g.55241707G>T")

## customize fields
getVariant("chr7:g.55241707G>T", 
          fields=c("dbnsfp.cadd.phred", "dbnsfp.polyphen2"), 
          return.as="text")
}


