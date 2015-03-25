### Various functions for reading VCF columns and creating HGVS IDs
library(VariantAnnotation)

getVcf <- function(file.path){
  snpVcf <- read.csv(file.path, stringsAsFactors=FALSE, header=F, sep='\t', comment.char="#")
  names(snpVcf) <- c("CHROM", "POS", "rsID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
  snpVcf
}

readSnps <- function(vcf){
  subs <- subset(vcf, nchar(REF) == nchar(ALT))
  snps <- subset(subs, nchar(REF) == 1)
  if(length(snps) > 0){
    hgvs <- paste(vcf$CHROM, ":g.", vcf$POS, vcf$REF, ">", vcf$ALT, sep="")
  }
  else {
    hgvs <- NULL
  }
  hgvs
}

readDels <- function(vcf){
  dels <- subset(vcf, nchar(REF) > nchar(ALT))
  if(length(dels) > 0){
    end <- dels$POS + (nchar(dels$REF) - 1)
    hgvs <- paste(dels$CHROM, ":g.", dels$POS,
                  "_", end, "del", sep="")
  }
  else {
    hgvs <- NULL
  }
  hgvs
}

readIns <- function(vcf){
  subs <- subset(vcf, nchar(REF) < nchar(ALT))
  if(length(subs) > 0){
    alt <- subs$ALT
    ins <- unlist(lapply(alt, function(i) substring(i, 2, nchar(i))))
    end <- subs$POS + nchar(ins)
    hgvs <- paste(subs$CHROM, ":g.", subs$POS,
                  "_", end, "ins", ins, sep="")
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
  hgvs <- list(snps=snps, deletions=dels, insertions=ins)#, indels=indels)
  hgvs
}

#df.list <- lapply(list(hgvs[1:10000],hgvs[10001:20000],hgvs[20001:30000], hgvs[30001:length(hgvs)]), getVariants)
#dfs <- lapply(df.list, as.data.frame)
#snps <- do.call(rbind.fill, dfs)

