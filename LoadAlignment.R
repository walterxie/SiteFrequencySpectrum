# Install b2xml
#library("devtools")
#devtools::install_github("walterxie/b2xml")
#library("b2xml")

setwd("~/WorkSpace/SiteFrequencySpectrum")
require(tidyverse)
require(xml2)

# load BEAST2 xml
xml.path = "~/WorkSpace/sc-analyses/xmls/CRC09_coal_simple.xml"
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
cat("Processed", length(seqs), "taxa, genotype = ", comma+1)
seqs[[1]][1:3]

nsite = length(seqs[[1]])
ntaxa = length(seqs)
ngt = comma+1

for (n in 1:nsite) {
  # per site or sample?
  # gtmat = matrix(ncol=ngt, nrow=nsite)
  gtmat = matrix(ncol=ngt, nrow=ntaxa)
  
  for (t in 1:ntaxa) {
    site <- seqs[[t]][n]
    # GT is split by , for a site 
    gt <- as.numeric(str_split(site, ",", simplify = T))
    
    # gtmat[n, ] <- gt
    gtmat[t, ] <- gt
  }
  
  # allele frequency 
  freq <- colSums(gtmat)
  
}


dat <- data.frame(gt=1:ngt, af=freq)
ggplot(dat, aes(x=factor(gt), y=af)) + 
  geom_bar(stat = "identity") +
  xlab("Genotype") + ylab("Allele frequency") +
  theme_bw()




