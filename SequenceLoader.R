

xmlToSequences <- function(xml.path) {
  require(tidyverse)
  require(xml2)
  
  b2xml = xml2::read_xml(xml.path)
  seqxml <- xml_find_all(b2xml, xpath="//sequence")
  
  comma = 0
  seqs <- list()
  for (s in seqxml) {
    # <sequence taxon="TI-TNS10" uncertain="true">
    taxon = xml_attrs(s)["taxon"]
    # rm new line char
    seqStr = gsub("[\r\n]", "", xml_text(s))
    # site is split by ;
    sites <- str_split(seqStr, ";", simplify = T)
    # rm last empty element because sequence has ; in the end
    sites <- sites[-length(sites)]
    
    cat("taxon ", taxon, " has", length(sites), "sites\n")
    
    if (comma < 1) {
      comma = str_count(sites, ",")[1]
    } 
    stopifnot(all(str_count(sites, ",")  == comma))
    
    seqs[[taxon]] <- sites
  }
  # gt = number of comma + 1
  cat("Processed", length(seqs), "taxa,", length(seqs[[1]]), "sites, which have", comma+1, "genotypes.")
  return(seqs)
}

# return a probability matrix for one site, rows are taxa, columns are GT 
getOneSiteProbUnphased <- function(siteID = 1, seqs=list()) {
  site <- sapply(seqs, `[`, siteID)
  # 1 site
  gtvec <- as.numeric(str_split(site, ",", simplify = T))
  # rows are taxa, columns are GT
  gtm <- matrix(gtvec, nrow = length(site)) 
  rownames(gtm) = names(site)    
  
  # convert to probability
  probPhased = prop.table(gtm, 1)
  # sapply(GT16, function(x){which(GT16unphased %in% x)})
  
  # convert to unphased
  # AA
  probUnphased = matrix(probPhased[,1], nrow=nrow(probPhased), ncol=1)
  # AC
  probUnphased = cbind(probUnphased, rowSums(probPhased[,c(2,5)]))
  probUnphased = cbind(probUnphased, rowSums(probPhased[,c(3,9)]))
  probUnphased = cbind(probUnphased, rowSums(probPhased[,c(4,13)]))
  # CC
  probUnphased = cbind(probUnphased, probPhased[,6])
  probUnphased = cbind(probUnphased, rowSums(probPhased[,c(7,10)]))
  probUnphased = cbind(probUnphased, rowSums(probPhased[,c(8,14)]))
  # GG
  probUnphased = cbind(probUnphased, probPhased[,11])
  probUnphased = cbind(probUnphased, rowSums(probPhased[,c(12,15)]))
  # TT
  probUnphased = cbind(probUnphased, probPhased[,16])
  
  return(probUnphased)
}



### in dev

# use sapply(seqs, `[`, 2) to get a site for all sample
likelihoodToProbPerSitePerSample <- function(site="", method=c("max","prob"),
                                             GT16=c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"),
                                             GT16unphased=c("AA","AC","AG","AT","AC","CC","CG","CT","AG","CG","GG","GT","AT","CT","GT","TT")) {
  # GT is split by , e.g. "0.1,0,2.5119e-01,0,0,0,0,0,2.5119e-01,0,7.9433e-04,0,0,0,0,0" 
  gtvec <- as.numeric(str_split(site, ",", simplify = T))
  # missing data or the character "?", e.g. "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1" 
  if ( all(gtvec== 1) ) {
    return("?")
  }
  ## evaluate choices
  method <- match.arg(method)
  
  if (method == "max") {
    # determine GT with highest prob 
    gt <- which(gtvec==max(gtvec))
    
    # combine unphased GT
    if (length(gt) > 1) {
      gt = which(GT16 %in% unique(GT16unphased[gt]))
    }
    
  } else {
    
    
    #TODO
    
  }
  
  stopifnot( !(is.null(gt) | length(gt) > 3 | length(gt) < 1) )
  # 1 or 2
  return(paste0(GT16[gt], collapse = "|"))
}



