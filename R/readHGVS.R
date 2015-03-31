### Various functions for reading VCF columns and creating HGVS IDs
library(S4Vectors)

getVcf <- function(file.path){
  Vcf <- read.csv(file.path, stringsAsFactors=FALSE, header=F, sep='\t', comment.char="#")
  names(Vcf) <- c("CHROM", "POS", "rsID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
  Vcf
}

readSnps <- function(vcf){
  subs <- subset(vcf, nchar(REF) == nchar(ALT))
  snps <- subset(subs, nchar(REF) == 1)
  if(nrow(snps) > 0){
    hgvs <- data.frame("query"=paste(snps$CHROM, ":g.", snps$POS, snps$REF, ">", snps$ALT, sep=""),
                       "type"=rep("snp", nrow(snps)), 
                       "pos"=paste(snps$CHROM, ":", snps$POS, "-", snps$POS, sep=""))
  }
  else {
    hgvs <- NULL
  }
  hgvs
}

readDels <- function(vcf){
  dels <- subset(vcf, nchar(REF) > nchar(ALT))
  if(nrow(dels) > 0){
    end <- dels$POS + (nchar(dels$REF) - 1)
    hgvs <- data.frame("query"=paste(dels$CHROM, ":g.", dels$POS,
                  "_", end, "del", sep=""), 
                  "type"=rep("deletion", nrow(dels)), 
                  "pos"=paste(dels$CHROM, ":", dels$POS, "-", dels$POS, sep=""))
  }
  else {
    hgvs <- NULL
  }
  hgvs
}

readIns <- function(vcf){
  subs <- subset(vcf, nchar(REF) < nchar(ALT))
  if(nrow(subs) > 0){
    alt <- subs$ALT
    ins <- unlist(lapply(alt, function(i) substring(i, 2, nchar(i))))
    end <- subs$POS + 1
    hgvs <- data.frame("query"=paste("id"=subs$CHROM, ":g.", subs$POS,
                  "_", end, "ins", ins, sep=""), 
                  "type"=rep("insertion", nrow(subs)), 
                  "pos"=paste(subs$CHROM, ":", subs$POS, "-", subs$POS, sep=""))
  }
  else {
    hgvs <- NULL
  }
  hgvs
}

# readIndels <- function(vcf){
#   subs <- subset(vcf, nchar(REF) > nchar(ALT))
#   indels <- subset(subs, nchar(ALT)) > 1)
#   if(length(indels) > 0){
#     hgvs <- paste(indels$CHROM, ":g.", indels$POS, 
#                   "_", indels$POS, indels$ALT, sep="")
#   }
#   else{ 
#     hgvs <- NULL}
#   hgvs
# }

getHgvs <- function(vcf){
  snps <- readSnps(vcf)
  dels <- readDels(vcf)
  ins <- readIns(vcf)
  #indels <- readIndels(vcf)
  hgvs <- do.call(rbind, list(snps, dels, ins))
  hgvs
}

.pasteChr <- function(hgvs.id){
  if (grepl("chr", hgvs.id)){
    return(hgvs.id)
  }
  else{
    hgvs <- paste("chr", hgvs.id, sep="")
    return(hgvs)
  }
}


