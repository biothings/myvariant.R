#library(myvariant)
#library(mygene)
library(magrittr)
library(S4Vectors)
library(VariantAnnotation)
library(plyr)

#functions
# return data.frame of Snp annotations
vcfSnpFilter <- function(vcf.file.input){
  vcf <- readVcf(vcf.file.input, genome="hg19")
  hgvs <- formatHgvs(vcf, "snp")
  mv.snps <- getVariants(hgvs)
  return(mv.snps)
}

# return data.frame of indel annotations
vcfIndel <- function(vcf.file.input){
  vcf <- getVcf(vcf.file.input)
  ## select high quality variants
  pass <- subset(vcf, FILTER == "PASS")
  ## get HGVS IDs and variant annotations
  hgvs <- getHgvs(pass)
  hgvs <- getHgvs(vcf)
  indels <- subset(hgvs, type!="snp")
  
  ## indels:
  if (nrow(indels) > 0){
    print(nrow(indels))
    indels$GENE <- sapply(indels$pos, function(i) query(i, fields="entrezgene", species="human")$hits$entrezgene)
    mv.indels <- getVariants(indels$query)
    indels <- merge(mv.indels, indels, by=intersect(names(indels), names(mv.indels)))
  }
  indels
}


replaceWith0 <- function(df){
  d <- data.frame(df)
  d[is.na(d)] <- 0
  DataFrame(d)
}

## apply filters to dataframe
filterDf <- function(df){
  df <- subset(df, cadd.consequence %in% c("NON_SYNONYMOUS", "STOP_GAINED", "STOP_LOST", "CANONICAL_SPLICE", "SPLICE_SITE"))
  df <- subset(df, exac.af < 0.01)
  #df <- subset(df, dbsnp.dbsnp_build > 128 )
  df <- subset(df, sapply(dbnsfp.1000gp1.af, function(i) i < 0.01 ))
  df
}

## get specific gene's rows of data.frame
geneInDf <- function(df, gene){
  gene.df <- subset(df, sapply(df$dbnsfp.genename, function(i) gene %in% i))
  gene.df
}

max.of.df <- function(df){
    #return(data.frame(subset(df, cadd.phred == max(df$cadd.phred))))
    return(data.frame(subset(df, cadd.phred == median(df$cadd.phred))))
}

## create unique dataframes by gene, check pathogenicity
geneDf <- function(vars.list, gene){
  patho <- lapply(vars.list, function(i) geneInDf(i, gene))
  patho <- patho[sapply(patho, function(i) nrow(i) > 0)]
  common <- do.call(rbind.fill, lapply(patho, max.of.df))
  df <- list("gene"=as.character(gene),
             "cadd.phred"=mean(common$cadd.phred))
  df
}

variantInDf <- function(df, variant){
      variant.df <- subset(df, sapply(df$Variant, function(i) variant %in% i))
      variant.df
}

variantDf <- function(vars.list, variant){
  var <- lapply(vars.list, function(i) variantInDf(i, variant))
  var
}

rankByCaddScore <- function(gene.list, df.list){
  y <- do.call(rbind, lapply(gene.list, function(i) geneDf(df.list, i)))
  df <- data.frame(gene=unlist(y[,1]), cadd.phred=unlist(y[,2]))
  ranked <- arrange(df, -cadd.phred)
  ranked
}
collapse <- function(...) {
  paste(unlist(list(...)), sep=",", collapse=",")
}
