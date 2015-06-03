### Various functions for reading VCF columns and creating HGVS IDs
library(S4Vectors)
library(plyr)
library(magrittr)

getHgvs <- function(chrom, pos, ref, alt, mutant.type=FALSE){
  if(nchar(ref) == nchar(alt) & nchar(ref) == 1){
    ## snp
    hgvs <- paste(chrom, ":g.", pos, ref, ">", alt, sep="")
    if(mutant.type){var.type <- "snp"}
  }
  else if(nchar(ref) > 1 & nchar(alt) == 1){
    ## deletion
    if(substr(ref, 1, 1) == alt){
      start <- as.integer(pos) + 1
      end <- as.integer(pos) + nchar(ref) - 1
      hgvs <- paste(chrom, ":g.", start, "_", end, "del", sep="")
      if(mutant.type){var.type <- "deletion"}
    }
    else{
      end <- as.integer(pos) + nchar(ref) - 1
      hgvs <- paste(chrom, ":g.", pos, "_", end, "delins", alt, sep="")
      if(mutant.type){var.type <- "delins"}
    }
  }
  else if(nchar(ref) == 1 & nchar(alt) > 1){
    ## insertion
    if(substr(alt, 1, 1) == ref){
      hgvs <- paste(chrom, ":g.", pos, "_", (as.integer(pos) + 1), "ins", sep="")
      ins_seq <- substr(alt, 2, nchar(alt))
      hgvs <- paste(hgvs, ins_seq, sep="")
      if(mutant.type){var.type <- "insertion"}
    }
    else{
      hgvs <- paste(chrom, ":g.", pos, "delins", alt, sep="")
      if(mutant.type){var.type <- "delins"}
    }
  }
  else if(nchar(ref) > 1 & nchar(alt) > 1){
    end <- as.integer(pos) + nchar(alt) - 1
    hgvs <- paste(chrom, ":g.", pos, "_", end, "delins", alt, sep="")
    if(mutant.type){var.type <- "delins"}
  }
  else{stop("Cannot convert pos, chrom, ref, alt into HGVS id")}
  if(!grepl("chr", hgvs)){
    hgvs <- paste("chr", hgvs, sep="")
  }
  if(mutant.type){
    c(hgvs, var.type)
  }
  else{hgvs}
}






getVcf <- function(file.path){
  stopifnot(grepl(".vcf", file.path))
  Vcf <- read.csv(file.path, stringsAsFactors=FALSE, header=FALSE, sep='\t', comment.char="#")
  names(Vcf) <- c("CHROM", "POS", "rsID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")[1:ncol(Vcf)]
  Vcf <- normalize.vcf(Vcf)
  if(!grepl("chr", Vcf$CHROM)){
    Vcf$CHROM <- paste("chr", Vcf$CHROM, sep="")
  }
  Vcf
#   Vcf <- DataFrame(transform(Vcf, colsplit(V10, split = "\\:", names = c('GT', 'AD', 'DP', 'GQ', 'PL')), stringsAsFactors=F))
#   Vcf$AD <- .factor2List(Vcf$AD)
#   Vcf$PL <- .factor2List(Vcf$PL)
}


getAll <- function(vcf.df){
  stopifnot(all(c('CHROM', 'POS', 'rsID', 'REF', 'ALT') %in% colnames(vcf.df)))
  snps <- getSnps(vcf.df)
  dels <- getDels(vcf.df)
  ins <- getIns(vcf.df)
  hgvs <- do.call(plyr::rbind.fill, list(snps, dels, ins))
  hgvs
}

getSnps <- function(vcf.df){
  stopifnot(all(c('CHROM', 'POS', 'rsID', 'REF', 'ALT') %in% colnames(vcf.df)))
  vcf.df %>%
    subset(!grepl(",", ALT) & nchar(REF) == nchar(ALT) &
             nchar(REF) == 1) %>%
    transform(query=paste(CHROM, ":g.", POS, REF, ">", ALT, sep=""),
              type="snp", 
              pos=paste(.trim(CHROM), ":", .trim(POS), "-", .trim(POS), sep=""))
}

getDels <- function(vcf.df){
  stopifnot(all(c('CHROM', 'POS', 'rsID', 'REF', 'ALT') %in% colnames(vcf.df)))
  vcf.df %>%
    subset(!grepl(",", ALT) & nchar(REF) > nchar(ALT) &
             substring(REF, 1, 1) == ALT) %>%
    transform(query=paste(CHROM, ":g.", POS,
                  "_", (POS + nchar(REF) - 1), "del", sep=""),
                  type="deletion",
                  pos=paste(.trim(CHROM), ":", .trim(POS), "-", .trim(POS), sep=""))
}

getIns <- function(vcf.df){
  stopifnot(all(c('CHROM', 'POS', 'rsID', 'REF', 'ALT') %in% colnames(vcf.df)))
  insertions <- subset(vcf.df, !grepl(",", ALT) & nchar(REF) < nchar(ALT) &
                  REF == substring(ALT, 1, 1))
  ins <- unlist(lapply(insertions$ALT, function(i) substring(i, 2, nchar(as.vector(i)))))
  end <- insertions$POS + 1
  hgvs <- cbind(insertions, data.frame(query=paste("id"=insertions$CHROM, ":g.", insertions$POS,
                  "_", end, "ins", ins, sep=""), 
                  type="insertion", 
                  pos=paste(.trim(insertions$CHROM), ":", .trim(insertions$POS), "-", .trim(insertions$POS), sep="")))
  hgvs
}

getIndels <- function(vcf.df){
  stopifnot(all(c('CHROM', 'POS', 'rsID', 'REF', 'ALT') %in% colnames(vcf.df)))
  vcf <- subset(vcf.df, !grepl(",", ALT))
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
                     pos=paste(.trim(delins$CHROM), ":", .trim(delins$POS), "-", .trim(delins$POS), sep="")))
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
