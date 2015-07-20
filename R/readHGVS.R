### Various functions for reading VCF columns and creating HGVS IDs
library(S4Vectors)
library(plyr)
library(magrittr)
library(VariantAnnotation)
library(GenomeInfoDb)

formatSingleHgvs <- function(chrom, pos, ref, alt, mutant_type=FALSE){
  if(nchar(ref) == nchar(alt) & nchar(ref) == 1){
    ## snp
    hgvs <- paste(chrom, ":g.", pos, ref, ">", alt, sep="")
    if(mutant_type){var.type <- "snp"}
  }
  else if(nchar(ref) > 1 & nchar(alt) == 1){
    ## deletion
    if(substr(ref, 1, 1) == alt){
      start <- as.integer(pos) + 1
      end <- as.integer(pos) + nchar(ref) - 1
      hgvs <- paste(chrom, ":g.", start, "_", end, "del", sep="")
      if(mutant_type){var.type <- "deletion"}
    }
    else{
      end <- as.integer(pos) + nchar(ref) - 1
      hgvs <- paste(chrom, ":g.", pos, "_", end, "delins", alt, sep="")
      if(mutant_type){var.type <- "delins"}
    }
  }
  else if(nchar(ref) == 1 & nchar(alt) > 1){
    ## insertion
    if(substr(alt, 1, 1) == ref){
      hgvs <- paste(chrom, ":g.", pos, "_", (as.integer(pos) + 1), "ins", sep="")
      ins_seq <- substr(alt, 2, nchar(alt))
      hgvs <- paste(hgvs, ins_seq, sep="")
      if(mutant_type){var.type <- "insertion"}
    }
    else{
      hgvs <- paste(chrom, ":g.", pos, "delins", alt, sep="")
      if(mutant_type){var.type <- "delins"}
    }
  }
  else if(nchar(ref) > 1 & nchar(alt) > 1){
    end <- as.integer(pos) + nchar(alt) - 1
    hgvs <- paste(chrom, ":g.", pos, "_", end, "delins", alt, sep="")
    if(mutant_type){var.type <- "delins"}
  }
  else{stop("Cannot convert pos, chrom, ref, alt into HGVS id")}
  if(!grepl("chr", hgvs)){
    hgvs <- paste("chr", hgvs, sep="")
  }
  if(mutant_type){
    c(hgvs, var.type)
  }
  else{hgvs}
}

formatHgvs <- function(vcf, variant_type=c("snp", "insertion", "deletion")){
  seqlevelsStyle(vcf) <- "UCSC"
  if ("snp" %in% variant_type){
    snps <- .getSnps(vcf)
  }
  else{ snps <- NULL}
  if ("insertion" %in% variant_type){
    ins <- .getIns(vcf)
  }
  else{ins <- NULL}
  if ("deletion" %in% variant_type){
    del <- .getDels(vcf)
  }
  else{del <- NULL}
  hgvs <- c(snps, ins, del)
  hgvs
}

.getSnps <- function(vcf){
  snp <- rowRanges(vcf)[isSNV(vcf)]
  if (length(snp) > 0){
  hgvs <- paste(seqnames(snp), ":g.", start(snp), as.character(snp$REF), ">", 
                as.character(unlist(snp$ALT)), sep="")
  }
  else{hgvs <- NULL}
  hgvs
}

.getDels <- function(vcf){
  del <- rowRanges(vcf)[isDeletion(vcf)]
  if (length(del) > 0){
  hgvs <- paste(seqnames(del), ":g.", start(del),
                  "_", end(del), "del", sep="")
  }
  else {hgvs <- NULL}
  hgvs             
}

.getIns <- function(vcf){
  ins <- rowRanges(vcf)[isInsertion(vcf)]
  if (length(ins) > 0) {
  hgvs <- paste(seqnames(ins), ":g.", start(ins),
                  "_", end(seq), "ins", alt(ins), sep="")
  }
  else {hgvs <- NULL}       
  hgvs
}

.getIndels <- function(vcf){
  vcf <- rowRanges(vcf)[isDelins(vcf)]
  ## case 1, nchar(ALT) == 1
  dels <- subset(vcf, nchar(REF) > 1 & nchar(ALT) == 1)
  hgvs.1 <- NULL
  if(nrow(dels) > 0){
  hgvs.1 <- paste(dels$CHROM, ":g.", dels$POS, 
                  "_", (dels$POS + nchar(dels$REF) - 1), "delins", dels$ALT, sep="")
  }
  ## case 2, nchar(REF) == 1
  ins <- subset(vcf, nchar(REF) == 1 & nchar(ALT) > 1)
  hgvs.2 <- NULL
  if(nrow(ins) > 0){
  hgvs.2 <- paste(ins$CHROM, ":g.", ins$POS,
                  "delins", ins$ALT)
  }
  ## case 3, 
  indel <- subset(vcf, nchar(REF) > 1 & nchar(ALT) > 1)
  hgvs.3 <- NULL
  if(nrow(indel) > 0){
  hgvs.3 <- paste(indel$CHROM, ":g.", indel$POS,
                  "_", (indel$POS + nchar(indel$ALT) - 1),
                  "delins", indel$ALT)
  }
  delins <- do.call(rbind, list(dels, ins, indel))
  if(nrow(delins) > 0){
  hgvs <- data.frame(query=c(hgvs.1, hgvs.2, hgvs.3),
                     type=rep("indel", nrow(dels) + nrow(ins) + nrow(indel),
                     pos=paste(.trim(delins$CHROM), ":", .trim(delins$POS), "-", 
                               .trim(delins$POS), sep="")))
  }
  else{hgvs <- NULL}
  hgvs
}

## normalizes rows where ALT == "GA,G" (multiple ALT values)
normalize.vcf <- function(vcf.df){
  if (nrow(vcf.df) == 0)
    return(vcf.df)
  split.alt <- strsplit(vcf.df$ALT, ",")
  vcf <- vcf.df[rep(seq(nrow(vcf.df)), elementLengths(split.alt)),]
  vcf$ALT <- unlist(split.alt)
  vcf
}


.trim <- function(x){
  gsub("^\\s+|\\s+$", "", x)
}
