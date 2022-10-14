# convert GT16 data in Phylonco XML into allele frequency

setwd("~/WorkSpace/SiteFrequencySpectrum")
require(tidyverse)
require(xml2)

### process BEAST2 xml
xml.path = "~/WorkSpace/sc-analyses/xmls/CRC09_coal_simple.xml"

seqs <- xmlToSequences(xml.path)
seqs[[1]][1:3]

nsite = length(seqs[[1]])
ntaxa = length(seqs)


### call SNPs

GT16 = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
# unphased should always collapse the heterozygotes to the alphabetic version (e.g CG|GC => CG)
GT16unphased = 
       c("AA","AC","AG","AT","AC","CC","CG","CT","AG","CG","GG","GT","AT","CT","GT","TT")

# last row must be healthy
healthyID = ntaxa
stopifnot( "healthy" == tolower(names(seqs)[healthyID]) )

gtmat = matrix(nrow=ntaxa, ncol=nsite)
for (n in 1:nsite) {
  # per site or sample?
  # gtmat = matrix(ncol=ngt, nrow=nsite)
  
  for (t in 1:ntaxa) {
    site <- seqs[[t]][n]
    # GT is split by , for a site 
    gtvec <- as.numeric(str_split(site, ",", simplify = T))
    # missing data or the character "?"
    if ( all(gtvec== 1) ) {
      # last row is healthy
      stopifnot( t != healthyID )
      gtmat[t, n] <- "?"
      next
    }
    # determine GT with highest prob 
    gt <- which(gtvec==max(gtvec))
    
    # combine unphased GT
    if (length(gt) > 1) {
      gt = which(GT16 %in% unique(GT16unphased[gt]))
    }
    
    stopifnot( !(is.null(gt) | length(gt) > 3 | length(gt) < 1) )
    # 1 or 2
    gtmat[t, n] <- paste0(GT16[gt], collapse = "|")
    
    # validate healthy
    if (t == healthyID) {
      # must be certain
      stopifnot( length(gt) == 1 )
    }
  }
}
# add taxa to 1st col
gtdf <- gtmat %>% as_tibble() %>% add_column(taxa = names(seqs), .before = 1)
print(gtdf, n =ntaxa)
gtdf %>% count(V1, sort = TRUE)
stringr::str_count(gtdf[[2]])

write_tsv(gtdf, file.path("data", "CRC09-SNPs-unphased.tsv"))

### minor allele frequency 

# rm healthy
snps <- gtdf %>% filter(row_number()!=healthyID)
# 25 × 24,386
print(snps[,1], n=26)
unique(unlist(apply(snps[,-1], 2, unique)))

### method 2 : sum the minimum counts
# e.g. CC=14, ?=8, CG=2, GG=2, then MA=2+2

# by column
countMinorAlleleSumMin <- function(snpv=c()) {
  require(tidyverse)
  af = snpv %>% tibble::as_tibble() %>% dplyr::count(value, sort = TRUE)
  # all discovered allele counts
  minaf = af %>% filter(value != "?") 
  # after rm "?"
  if (nrow(minaf) <= 1 || min(minaf$n)==max(minaf$n)) {
    # uncertain minor alleles
    # e.g. AA=8,AG=8,GG=8, or only 1 or no GT left
    return(c(paste0(minaf[["value"]], collapse = "/"), 0))
  }
  # take all min(n)
  minaf = minaf %>% filter(n == min(n))
  # sum all min(n)
  return( c(paste0(minaf[["value"]], collapse = "/"), 
            sum(minaf[["n"]])) )
}
names(gtdf[,-1])[1:10]
# countMinorAlleleSumMin(gtdf[["V1"]])
mac = lapply(gtdf[,-1], countMinorAlleleSumMin)
mac[1:3]
stopifnot(length(mac)==ncol(gtdf)-1)

mac.df <- as_tibble(mac) %>% as.matrix %>% t %>% as_tibble
names(mac.df) = c("minorallele","count")

write_tsv(mac.df, file.path("data", "CRC09-minor-allele-count-with-healthy.tsv"))

### re-load saved data 
require(tidyverse)
require(ggplot2)

setwd("~/WorkSpace/SiteFrequencySpectrum")
mac.df = read_tsv(file.path("data", "CRC09-minor-allele-count-with-healthy.tsv"))

# summary 
unique(mac.df[[1]])
mac.df[[2]] = as.numeric(mac.df[[2]])
# counts
sort(unique(mac.df[[2]]))
# mac.df %>% filter(count == 16)
# which(mac.df$count == 16)
# snps[["V807"]] %>% tibble::as_tibble() %>% dplyr::count(value, sort = TRUE)
# gtdf[[807+1]][healthyID] # id+1
# mac.df[807,]
# problems, site: 5, 807, 2578

# plot 
ntaxa=26
# exclude healthy
print(ntaxa)

mac.df[["count"]] = as.numeric(mac.df[["count"]])
p3 <- ggplot(mac.df %>% filter(count > 0), aes(count)) + 
  geom_histogram(stat = 'count') + 
  scale_x_continuous(breaks = sort(unique(mac.df[["count"]]))) +
  xlab(paste("Minor allele count in a site for", (ntaxa), "samples")) + 
  theme_bw()

ggsave(file.path("figures", "minor-allele-count-spectrum-with-healthy.pdf"), p3)

### method 3 : rm ? and 1 of maximun count(s), and then sum the rest 

countMinorAlleleSumTheRest <- function(snpv=c()) {
  require(tidyverse)
  af = snpv %>% tibble::as_tibble() %>% dplyr::count(value, sort = TRUE)
  # all discovered allele counts
  nounk = af %>% filter(value != "?") 
  # after rm "?", 
  if (nrow(nounk) <= 1) {
    # uncertain minor alleles
    # e.g. AA=8,AG=8,GG=8, or only 1 or no GT left
    return(c("", 0))
  }
  
  # rm max
  maxid = which(nounk$n == max(nounk$n))
  # if max==min
  if (length(maxid) >= 1) {
    # try to rm 1 of max
    minaf = nounk[-maxid[1],]
    # sum all min(n)
    return( c(paste0(minaf[["value"]], collapse = "/"), 
              sum(minaf[["n"]])) )
  } else {
    return(c("", 0))
  }
}
# with healthy
# countMinorAlleleSumTheRest(gtdf[["V1"]])
mac = lapply(gtdf[,-1], countMinorAlleleSumTheRest)
mac[1:3]
stopifnot(length(mac)==ncol(gtdf)-1)

mac.df <- as_tibble(mac) %>% as.matrix %>% t %>% as_tibble
names(mac.df) = c("minorallele","count")

write_tsv(mac.df, file.path("data", "CRC09-minor-allele-sum-rest-with-healthy.tsv"))

unique(mac.df[[1]])
sort(unique(mac.df[[2]]))

mac.df[["count"]] = as.numeric(mac.df[["count"]])
p4 <- ggplot(mac.df %>% filter(count > 0), aes(count)) + 
  geom_histogram(stat = 'count') + 
  scale_x_continuous(breaks = sort(unique(mac.df[["count"]]))) +
  xlab(paste("Minor allele count in a site for", (ntaxa), "samples")) + 
  theme_bw()
p4

ggsave(file.path("figures", "minor-allele-sum-rest-with-healthy.pdf"), p4)



### method 1 : too slow

# rm healthy
gtm <- gtdf %>% filter(row_number()!=healthyID)
print(gtm[,1], n=26)
unique(unlist(apply(gtm[,-1], 2, unique)))
       
afdf = NULL
# 24,385 sites
for (n in 1:nsite) {
  colnm = paste0("V", n) 
  af = gtm %>% count(!!as.symbol(colnm), sort = TRUE)
  # all discovered allele counts
  majoraf = af %>% filter(!!as.symbol(colnm) != "?" & n > min(n))
  minoraf = af %>% filter(n == min(n))
  if (majoraf %>% nrow == 0) {
    # e.g. V2152 CT == TT == 8
    # take one as minor, and rests are major, not affecting counts
    coun = minoraf %>% select(n) %>% unlist
    alle = minoraf %>% select(!!as.symbol(colnm)) %>% unlist
    mcount = coun[1]
    majcount = sum(coun[2:length(coun)])
    
    minor = alle[1]
    major = paste0(alle[2:length(coun)], collapse="|")
    
  } else {
    majcount = majoraf %>% select(n) %>% sum %>% as.numeric
    major = majoraf %>% select(!!as.symbol(colnm)) %>% unlist %>% paste0(collapse="|")
    
    # TODO drop rows having uncertainty e.g. 7th site CC|CT 1
    maf = minoraf %>% filter(!str_detect(!!as.symbol(colnm), '\\|'))
    if (nrow(maf) == 0) {
      maf = minoraf  
    } 
    if (nrow(maf) != nrow(minoraf)) {
      cat("Warning : drop the count below at site ", n, "\n")
      print(minoraf %>% filter(str_detect(!!as.symbol(colnm), '\\|')))
    }
    # sum count if multiple minor alleles
    mcount = maf %>% select(n) %>% sum %>% as.numeric 
    # minor allele
    minor = maf %>% select(!!as.symbol(colnm)) %>% unlist %>% paste0(collapse="|")
  }
  # ? counts
  unknown = af %>% filter(!!as.symbol(colnm) == "?") %>% select(n) %>% as.numeric
  
  if (n == 1)
    afdf = tibble(minorallele = minor, count = mcount, majorallele = major, count2 = majcount, unknown = unknown)
  else
    afdf = afdf %>% add_row(minorallele = minor, count = mcount, majorallele = major, count2 = majcount, unknown = unknown)
  
  if (n %% 10000 == 0) {
    cat("Site ", n, " adding ", nrow(afdf), " rows :\n")
    print(afdf %>% slice(n())) 
  }
}
### slow end
afdf
stopifnot(nrow(afdf) == nsite)

write_tsv(afdf, file.path("data", "CRC09-allele-frequency.tsv"))


### figures

require(tidyverse)
require(ggplot2)

setwd("~/WorkSpace/SiteFrequencySpectrum")
afdf = read_tsv(file.path("data", "CRC09-allele-frequency.tsv"))

ntaxa=26
# exclude healthy
print(ntaxa-1)
unknowndf = afdf %>% select(unknown) %>% mutate(freq = unknown / (ntaxa-1)) 

p1 <- ggplot(unknowndf, aes(freq)) + 
  geom_histogram() + 
  xlab(paste("Unknown allele frequency in a site for", (ntaxa-1), "samples")) + 
  theme_bw()

ggsave(file.path("figures", "unknown-allele-frequency-spectrum.pdf"), p1)

# minor allele frequency
minorafdf = afdf %>% select(minorallele, count, unknown) %>% mutate(freq = count / (ntaxa-1-unknown)) 
minorafdf
unique(minorafdf %>% select(minorallele))

p2 <- ggplot(minorafdf, aes(freq)) + 
  geom_histogram() + 
  xlab(paste("Minor allele frequency in a site for", (ntaxa-1), "samples")) + 
  theme_bw()

ggsave(file.path("figures", "minor-allele-frequency-spectrum.pdf"), p2)

  

### in dev 
p <- ggplot(afdf, aes(x=factor(gt), y=af)) + 
  geom_bar(stat = "identity") +
  xlab("Genotype") + ylab("Allele frequency") +
  theme_bw()
