setwd("~/WorkSpace/SiteFrequencySpectrum")
source("SequenceLoader.R")


xml.path = "~/WorkSpace/sc-analyses/xmls/CRC09_coal_simple.xml"

seqs = xmlToSequences(xml.path)
seqs[[1]][1:5]

nsite = length(seqs[[1]])
ntaxa = length(seqs)

### likelihood => prob

# the number of data per site having prob >= 0.5
getNConfPerSite <- function(siteID = 1, seqs=list()) {
  site <- getOneSiteProbUnphased(siteID, seqs) 
  enoughPr = apply(site, 1, function(x){max(x) >= 0.5})
  # count TRUE
  return(sum(enoughPr))
}

nT = sapply(1:nsite, getNConfPerSite, seqs)

#apply(site, 1, function(x){which(x==max(x))})

p1 <- ggplot() + aes(nT)+ geom_histogram(binwidth=1, colour="black", fill="white") +
  scale_x_continuous(breaks = sort(unique(nT))) +
  xlab(paste("The spectrum of #samples (with GT Pr >= 0.5) in a site, total samples = ", ntaxa)) + 
  theme_bw()

ggsave(file.path("figures", "count-spectrum-sample-prob-0.5.pdf"), p1)


getMaxProbPerSite <- function(siteID = 1, seqs=list()) {
  site <- getOneSiteProbUnphased(siteID, seqs) 
  maxP = apply(site, 1, function(x){max(x)})
  return(maxP)
}

maxPr = lapply(1:nsite, getMaxProbPerSite, seqs)

frac = sapply(maxPr, function(x){ length(x[x>=0.5]) / length(x) })

p2 <- ggplot() + aes(frac)+ geom_histogram(colour="black", fill="white") +
  xlab(paste("The spectrum of fraction of #samples (with GT Pr >= 0.5) in a site, total samples = ", ntaxa)) + 
  theme_bw()
ggsave(file.path("figures", "fraction-spectrum-sample-prob-0.5.pdf"), p2)
