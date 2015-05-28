\name{metadata}
\docType{methods}
\alias{metadata}
\alias{metadata,MyVariant-method}

\title{
metadata
}
\description{
Get metadata for MyVariant.info services.
}
\usage{
metadata(x, ...)
}
\arguments{
  \item{x}{
MyVariant object
}
  \item{...}{
MyVariant object slot parameters
}}

\value{returns the metadata including available databases and number of documents.}

\references{
http://myvariant.info/v1/metadata
}

\examples{
## Get metadata 
myvariant<-MyVariant()
metadata(myvariant)

}
