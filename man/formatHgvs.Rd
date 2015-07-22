\name{formatHgvs}
\docType{methods}
\alias{formatHgvs}

\title{
Get all HGVS IDs from Vcf object.
}
\description{
Read in a Vcf object created by \code{\link{readVcf}} to extract all HGVS IDs for querying MyVariant.info.
}
\usage{
formatHgvs(vcf, variant_type = c("snp", "insertion", "deletion"))
}
\arguments{
  \item{vcf}{
Vcf object created by \code{\link{readVcf}} or VRanges object.
}
  \item{variant_type}{
Type of variant HGVS IDs to retrieve from Vcf object. Default c("snp", "insertion", "deletion")
}
}
\value{vector}

\references{
https://myvariant.info
http://www.hgvs.org/mutnomen/recs-DNA.html
}

\seealso{
\code{\link{formatSingleHgvs}}
}
\examples{
## return HGVS IDs for all snps in a Vcf
file.path <- system.file("extdata", "dbsnp_mini.vcf", package="myvariant")
vcf <- readVcf(file.path, genome="hg19")
hgvs <- formatHgvs(vcf, variant_type="snp")
}


