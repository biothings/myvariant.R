\name{getHgvs}
\docType{methods}
\alias{getHgvs}

\title{
Get Hgvs HGVS ID from chromosome, position, reference and alternate alleles.
}
\description{
Create a single HGVS ID for a variant from chromosome, position, reference and alternate alleles.
}
\usage{
getHgvs(chrom, pos, ref, alt, mutant.type=FALSE)
}
\arguments{
  \item{chrom}{
Chromosome.
}
  \item{pos}{
Position of the variant on the reference genome (hg19).
}
  \item{ref}{
Reference allele.
}
  \item{alt}{
Alternate allele.
}
  \item{mutant.type}{
Logical indicating whether to return the type of mutation along with the HGVS ID.
}
}
\value{returns a string}

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
## return HGVS ID for a variant
getHgvs(1, 35367, "G", "A")
}


