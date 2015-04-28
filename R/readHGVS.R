### Various functions for reading VCF columns and creating HGVS IDs
library(S4Vectors)
library(plyr)
library(magrittr)

getVcf <- function(file.path){
  Vcf <- read.csv(file.path, stringsAsFactors=FALSE, header=F, sep='\t', comment.char="#")
  names(Vcf) <- c("CHROM", "POS", "rsID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
  Vcf
}


getAll <- function(vcf){
  #vcf <- subset(vcf, FILTER=="PASS")
  vcf <- .normalize.subs(vcf)
  snps <- readSnps(vcf)
  dels <- readDels(vcf)
  ins <- readIns(vcf)
  #indels <- readIndels(vcf)
  hgvs <- do.call(rbind.fill, list(snps, dels, ins))#, indels))
  #hgvs$query <- lapply(.pasteChr(hgvs$query))
  hgvs
}

getSnps <- function(vcf){
  vcf %>%
    subset(!grepl(",", ALT) & nchar(REF) == nchar(ALT) &
             nchar(REF) == 1) %>%
    transform(query=paste(CHROM, ":g.", POS, REF, ">", ALT, sep=""),
              type="snp", 
              pos=paste(.trim(CHROM), ":", .trim(POS), "-", .trim(POS), sep=""))
}

getDels <- function(vcf){
  vcf %>%
    subset(!grepl(",", ALT) & nchar(REF) > nchar(ALT) &
             substring(REF, 1, 1) == ALT) %>%
    transform(query=paste(CHROM, ":g.", POS,
                  "_", (POS + nchar(REF) - 1), "del", sep=""),
                  type="deletion",
                  pos=paste(.trim(CHROM), ":", .trim(POS), "-", .trim(POS), sep=""))
}

getIns <- function(vcf){
  insertions <- subset(vcf, !grepl(",", ALT) & nchar(REF) < nchar(ALT) &
                  REF == substring(ALT, 1, 1))
  ins <- unlist(lapply(insertions$ALT, function(i) substring(i, 2, nchar(as.vector(i)))))
  end <- insertions$POS + 1
  hgvs <- data.frame(query=paste("id"=insertions$CHROM, ":g.", insertions$POS,
                  "_", end, "ins", ins, sep=""), 
                  type="insertion", 
                  pos=paste(.trim(insertions$CHROM), ":", .trim(insertions$POS), "-", .trim(insertions$POS), sep=""))
  hgvs
}

getIndels <- function(vcf){
  vcf <- subset(vcf, !grepl(",", ALT))
  ## case 1, nchar(ALT) == 1
  dels <- subset(vcf, nchar(REF) > 1 && nchar(ALT) == 1)
  hgvs.1 <- paste(dels$CHROM, ":g.", dels$POS, 
                  "_", (dels$POS + nchar(dels$REF) - 1), "delins", dels$ALT, sep="")
  
  ## case 2, nchar(REF) == 1
  ins <- subset(vcf, nchar(REF) == 1 & nchar(ALT) > 1)
  hgvs.2 <- paste(ins$CHROM, ":g.", ins$POS,
                  "delins", ins$ALT)
  
  ## case 3, 
  indel <- subset(vcf, nchar(REF) > 1 & nchar(ALT) > 1)
  hgvs.3 <- paste(indel$CHROM, ":g.", indel$POS,
                  "_", (indel$POS + nchar(indel$ALT) - 1),
                  "delins", indel$ALT)
  
  delins <- do.call(rbind, c(dels, ins, indel))
  hgvs <- data.frame(query=c(hgvs.1, hgvs.2, hgvs.3),
                     type=rep("indel", nrow(dels) + nrow(ins) + nrow(indel),
                     pos=paste(.trim(delins$CHROM), ":", .trim(delins$POS), "-", .trim(delins$POS), sep="")))
  hgvs
}

## normalizes rows where ALT == "GA,G" (multiple ALT values)
normalize.vcf <- function(vcf){
  if (nrow(vcf) == 0)
    return(vcf)
  split.alt <- strsplit(vcf$ALT, ",")
  vcf <- vcf[rep(seq(nrow(vcf)), elementLengths(split.alt)),]
  vcf$ALT <- unlist(split.alt)
  vcf
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

.trim <- function(x){
  gsub("^\\s+|\\s+$", "", x)
}
