### Various functions for reading VCF columns and creating HGVS IDs
library(VariantAnnotation)

readSnps <- function(vcf.object){
  subs <- subset(vcf.object, nchar(as.character(REF)) == nchar(as.character(unlist(ALT))))
  snps <- subset(subs, nchar(as.character(REF)) == 1)
  hgvs <- paste("chr", seqnames(snps), ":g.", start(snps), 
                as.character(snps$REF), ">", 
                as.character(unlist(snps$ALT)), sep="")
  hgvs
}

readDels <- function(vcf.object){
  dels <- subset(vcf.object, nchar(as.character(REF)) > nchar(as.character(unlist(ALT))))
  hgvs <- paste("chr", seqnames(dels), ":g.", `+`(start(dels), 1),
                "_", end(dels), "del", sep="")
  hgvs
}

readIns <- function(vcf.object){
  subs <- subset(vcf.object, nchar(as.character(REF)) < nchar(as.character(unlist(ALT))))
  alt <- as.character(unlist(subs$ALT))
  ins <- unlist(lapply(alt, function(i) substring(i, 2, nchar(i))))
  hgvs <- paste("chr", seqnames(subs), ":g.", `+`(start(subs), 1),
                "_", end(subs), "ins", ins, sep="")
  hgvs
}

getHgvs <- function(vcf.file){
  read <- readVcf(vcf.file, genome="hg19")
  vcf <- rowData(read)
  snps <- readSnps(vcf)
  dels <- readDels(vcf)
  ins <- readIns(vcf)
  hgvs <- c(snps, dels, ins)
  hgvs
}

#hgvs <- getHgvs("~/Documents/Su_Lab/myvariant.info/snp_analysis/miller.vcf")
#df.list <- lapply(list(hgvs[1:10000],hgvs[10001:20000],hgvs[20001:30000], hgvs[30001:length(hgvs)]), getVariants)
#dfs <- lapply(df.list, as.data.frame)
#snps <- do.call(rbind.fill, dfs)

