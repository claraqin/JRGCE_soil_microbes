# R code for Qin et al., Journal of Ecology, accepted Aug. 28, 2019

# Working directory must be set to source file's location for the 
# following code to read files properly.

RAREFACTION_ITS <- 28000
RAREFACTION_BAC <- 33500

# The script contains code for loading libraries, loading datasets,
# and calculating new metrics prior to formal statistical analysis.

###################
#  Load Packages  #
###################

library('phyloseq')
library('ggplot2')
library('vegan')
library('ecodist')
library('reshape2')
library('plyr')
library('dplyr')
library('leaps')
library('scales')
library('corrplot')
library('mvabund')
library('lme4')
library('lmerTest')
library('tidyr')
library('gdm')
library('MuMIn')
library('lemon')
theme_set(theme_bw())

##########################
# Load Plant & Soil Data #
##########################

# Load plant community data
Plants <- read.csv("../data/plant_biomass.csv")
Plants <- Plants[which(Plants$yr==2014),]
litter <- Plants[which(Plants$sp=="LITR"),]
Plants <- Plants[-which(Plants$sp=="LITR"),] # Remove litter
plants.w <- dcast(Plants, plt ~ sp, value.var="bio")
checksum <- colSums(plants.w)
plants.prune <- plants.w[,-c(which(checksum==0), 1)]
rownames(plants.prune) <- plants.w$plt
Plants.identified <- plants.prune[,!grepl("UN", names(plants.prune))]

plant.rich.all <- apply(plants.prune, 1, function(x) { sum(x!=0 & !grepl("UN", names(plants.prune)))})

# Load general sample data (includes soil data)
sampleData <- read.csv("../data/GCE_compiled_data.csv")
names(sampleData) <- tolower(names(sampleData))
for(col in 8:14) sampleData[,col]  <- as.factor(sampleData[,col])
for(col in c(15:16)) sampleData[,col] <- as.numeric(as.character(sampleData[,col]))
sampleData %>%
  dplyr::rename(lon = x,    # These are custom-scaled coordinates based on aerial imagery
                lat = y) ->
  sampleData

# Get Shannon diversity
shannon.div = function(spec.vec) {
  total.species = sum(spec.vec)
  normalized.vec0 = (spec.vec / total.species)
  normalized.vec = normalized.vec0[which(normalized.vec0!=0)]
  sum = 0
  for(p in normalized.vec) {
    sum = sum + p*log(p)
  }
  -sum
}

sampleData$plant.shannon <- apply(Plants.identified, 1, shannon.div)

# Remove extremely outlying perc.c and perc.n in sample GCE_7
sampleData$perc.c[which(sampleData$id==7)] <- NA
sampleData$perc.n[which(sampleData$id==7)] <- NA

# Correct for Carlo Erba batch effect by adjusting means of soil chem measurements
# First make sure that the CE batches are independent of (random with respect to)
# treatments, i.e. approx 50-50 ambient vs. elevated
output <- matrix(rep(NA,24), nrow=4)
rownames(output) <- seq(1,4)
colnames(output) <- c("n","co2","heat","precip","wildfire.2003","burn.2011")
for(i in 1:4) {
  output[i,1] <- mean(sampleData$n[which(sampleData$ce.batch==i)]==1)
  output[i,2] <- mean(sampleData$co2[which(sampleData$ce.batch==i)]==1)
  output[i,3] <- mean(sampleData$heat[which(sampleData$ce.batch==i)]==1)
  output[i,4] <- mean(sampleData$precip[which(sampleData$ce.batch==i)]==1)
  output[i,5] <- mean(sampleData$wildfire.2003[which(sampleData$ce.batch==i)]==1)
  output[i,6] <- mean(sampleData$burn.2011[which(sampleData$ce.batch==i)]==1)
}
output
# Approximately random
# Adjust means to match global mean
colrange <- c("perc.n","perc.c","n.mg","c.mg") # re-calculate c/n later
col3 <- which(names(sampleData)=="ce.batch")
correct_batch_means <- function(sample_data, colrange, col.batch) {
  global.means <- apply(sample_data[,colrange], 2, mean, na.rm=TRUE)
  for(i in 1:4) {
    ind <- which(sample_data[,col.batch]==i)
    batch.means <- apply(sample_data[ind,colrange], 2, mean, na.rm=TRUE)
    mean.diffs <- global.means - batch.means
    colrange.numeric <- match(colrange, names(sample_data))
    for(j in colrange.numeric) {
      sample_data[which(sample_data[,col.batch]==i),j] <- sample_data[which(sample_data[,col.batch]==i),j] + mean.diffs[j-min(colrange.numeric)+1]
    }
  }
  sample_data
}
sampleData <- correct_batch_means(sampleData, colrange, col3)

# C/N ratio
sampleData$c.n <- sampleData$perc.c / sampleData$perc.n


############################
#  Load ITS (Fungal) Data  #
############################

# Load data via phyloseq
ITS <- import_biom("../data/otu_table_wTax_ITS.biom",parseFunction=parse_taxonomy_greengenes)

# Match OTU table's sample IDs to sampleData's IDs
sampleDataITS <- import_qiime_sample_data("../data/mapping_file_ITS.txt")
sample_data(ITS) <- sampleDataITS
rm(sampleDataITS)

# Remove low-abundance taxa and poorly-classified taxa
ITS.prune <- prune_taxa(taxa_sums(ITS)>10, ITS)
kingdom <- as.vector(data.frame(tax_table(ITS))$Kingdom)
tax.keep <- kingdom=="Fungi"
tax.keep[is.na(tax.keep)] <- FALSE
ITS.prune <- prune_taxa(tax.keep,ITS)

# Histogram of sequencing depth
hist(sample_sums(ITS), col="orange", main="Sample sequencing depth\nITS")

# Which samples would be missing after rarefaction?
excluded.samples <- as.numeric(sub("F_GCE","",names(which(sample_sums(ITS.prune) < RAREFACTION_ITS))))
# Missing from which treatment combinations?
table(sampleData$treatment.group[match(sampleData$id, excluded.samples)])

# Rarefaction
ITS.rare <- rarefy_even_depth(physeq=ITS.prune,sample.size=RAREFACTION_ITS,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)
rm(ITS.prune)

###############################
#  Load 16S (Bacterial) Data  #
###############################

# Load data via phyloseq
BAC <- import_biom("../data/otu_table_wTax_16S.biom",parseFunction=parse_taxonomy_greengenes)

# Match OTU table's sample IDs to sampleData's IDs
sampleDataBAC <- import_qiime_sample_data("../data/mapping_file_16S.txt")
sample_data(BAC) <- sampleDataBAC
rm(sampleDataBAC)

# Drop the final column in the taxonomy table, "Rank1"
rm.ind <- which(colnames(tax_table(BAC)) == "Rank1")
tax_table(BAC) <- tax_table(BAC)[,-rm.ind]

# Also, remove sample GCE70 because it contained two missing values (perc.N and perc.C)
rm.ind <- which(rownames(sample_data(BAC))=="B_GCE70")
sample_data(BAC) <- sample_data(BAC)[-rm.ind,]

# Remove low-abundance taxa and poorly-classified taxa
BAC.prune <- prune_taxa(taxa_sums(BAC)>10, BAC)
kingdom <- as.vector(data.frame(tax_table(BAC))$Kingdom)
tax.keep <- kingdom=="Bacteria"
tax.keep[is.na(tax.keep)] <- FALSE
BAC.prune <- prune_taxa(tax.keep,BAC)

# Histogram of sequencing depth
hist(sample_sums(BAC), col="orange", main="Sample sequencing depth\nBAC")

# Which samples would be missing after rarefaction?
excluded.samples <- as.numeric(sub("B_GCE","",names(which(sample_sums(BAC.prune) < RAREFACTION_BAC))))
length(excluded.samples)
# Missing from which treatment combinations?
table(sampleData$treatment.group[match(sampleData$id, excluded.samples)])

# Rarefaction
BAC.rare <- rarefy_even_depth(physeq=BAC.prune,sample.size=RAREFACTION_BAC,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)
rm(BAC.prune)

#########################################
#### MERGE ALL DATA: sampleData.long ####
#########################################

rich.ITS <- cbind(estimate_richness(ITS.rare, split = TRUE, measures = c("Observed","Shannon")),
                  id = as(sample_data(ITS.rare),"data.frame")[,2])
names(rich.ITS) <- tolower(names(rich.ITS))
rich.ITS <- rich.ITS[order(rich.ITS$id),]

rich.BAC <- cbind(estimate_richness(BAC.rare, split = TRUE, measures = c("Observed","Shannon")),
                  id = as(sample_data(BAC.rare),"data.frame")[,2])
names(rich.BAC) <- tolower(names(rich.BAC))
rich.BAC <- rich.BAC[order(rich.BAC$id),]

# Get NMDS scores for plant community composition
# 3 NMDS axes; use first 2
set.seed(1)
Plants.ord <- vegan::metaMDS(Plants.identified, "bray", k=3)
stressplot(Plants.ord) # Non-metric fit 0.958, linear fit 0.744
Plants.scores <- scores(Plants.ord)
sampleData$NMDS.Plants.1 <- Plants.scores[match(sampleData$id, rownames(Plants.scores)), 1]
sampleData$NMDS.Plants.2 <- Plants.scores[match(sampleData$id, rownames(Plants.scores)), 2]

# Add plant richness to sampleData
sampleData$plant.rich <- plant.rich.all[match(sampleData$id, names(plant.rich.all))]

# Merge data
sampleData$plt <- as.factor(sampleData$plot)
merge(rich.ITS, rich.BAC, by="id", all=TRUE) %>% 
  dplyr::rename(ITS.observed = observed.x, ITS.shannon = shannon.x,
                BAC.observed = observed.y, BAC.shannon = shannon.y) %>%
  dplyr::select(id, ITS.observed, ITS.shannon, BAC.observed, BAC.shannon) %>%
  right_join(sampleData, by="id") %>%
  mutate(npp.std = (npp - mean(npp))/sd(npp)) %>%
  dplyr::select(id, n, co2, heat, precip, burn.2011, wildfire.2003, perc.n, c.n, water.content, ph, npp.std, plant.rich, plant.shannon, ITS.observed, ITS.shannon, BAC.observed, BAC.shannon, NMDS.Plants.1, NMDS.Plants.2, plt,
                contains("NMDS")) %>% # -> rich.microbial
  gather("variable","value",c("ITS.observed","ITS.shannon","BAC.observed","BAC.shannon")) ->
  sampleData.long
