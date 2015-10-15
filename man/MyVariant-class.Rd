\name{MyVariant-class}
\Rdversion{0.99}
\docType{class}
\alias{MyVariant-class}


\title{Class \code{"MyVariant"}}
\description{
R Client to access MyVariant.Info annotation services
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{MyVariant(base.url="http://myvariant.info/v1", delay=1,     step=1000, version=version, verbose=TRUE, debug=FALSE)}.

}
\section{Slots}{
  \describe{
    \item{\code{base.url}:}{"http://myvariant.info/v1". Object of class \code{"character"}}
    \item{\code{delay}:}{Sleep time between batch retrieval. Object of class \code{"numeric"}}
    \item{\code{step}:}{Batch limit. Object of class \code{"numeric"}}
    \item{\code{version}:}{httr package version. Object of class \code{"character"}}
    \item{\code{verbose}:}{Object of class \code{"logical"}}
    \item{\code{debug}:}{Object of class \code{"logical"}}
  }
}
\section{Methods}{
  \describe{
  \item{\code{getVariant(hgvsid, fields=NULL, 
        ..., return.as=c("records", "text"))}:}{ Return the variant object for the given hgvsid }
  \item{\code{getVariants(hgvsids, fields=NULL, ..., 
             return.as=c("DataFrame", "records", "text"))}:}{ Return the list of variant object for the given list of hgvsids. }
  \item{\code{queryVariant(q, fields=NULL, 
             ..., return.as=c("DataFrame", "records", "text"))}:}{ Return  the query result. }
  \item{\code{queryVariants(qterms, scopes=NULL, fields=NULL, 
          ..., return.as=c("DataFrame", "records", "text"), returnall=FALSE)}:}{ Return the batch query result. }
  }
}
\references{
https://github.com/SuLab/myvariant.info/wiki
}
\author{
Adam Mark
}

\examples{
showClass("MyVariant")
}
\keyword{classes}


