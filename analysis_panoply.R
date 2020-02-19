#
# analysis_panoply.R
#
# Example of analysis of the young/old dataset using panoplyCF
#
# 2020-02-18  WTR
#

library(wadeTools)
library(panoplyCF)

# edit these lines to reflect where you put stuff
#
#    data_base is where the original raw data are located. It's also where
#    you created the spreadsheet that contains the original demographics data
#    from FlowRepository, with the "valid" column added (in which you found
#    the bad instance, and set valid = 0 in that row).  So if you haven't done
#    that, do it now [hint - you can use excel, and save the result as csv].
#
#    gated_base is the directory that contains the gated FCS data
#
proj_base = "~/Data/Independent_Consulting/Penn/Matei/"
data_base = tight(proj_base, "data/young_old/")
gated_base = tight(proj_base, "results/young_old/gated_fcs/")

# use the spreadsheet to retrieve the data
tab = read.csv(tight(data_base, "Aging1_2_Demographics_validity.csv"))

# censor the bad instance.  You should have 135 rows after censoring
tab = tab[which(tab$valid == 1), ]

# make an aggregate of the data to train panoply
nper = 10000        # number of cell per instance
ninst = nrow(tab)   # number of instances in the dataset

set.seed(137)       # set the seed for reproducible sampling
fl = list()         # make a place to store the sample instances
for (i in 1:ninst) {
  cat("sampling", i, "of", ninst, "...")
  fn = tight(gated_base, tab$FCS.File[i])
  # sampling the flowFrames.  The easiest way is this:
  #     fl[[i]] = read.FCS(fn, which.lines = nper)
  # However, the following 2 lines are MUCH faster:
  tmp = read.FCS(fn)
  fl[[i]] = Subset(tmp, sampleFilter(size = nper))
  cat("done.\n")
}
# turn the list into a flowSet
fs = flowSet(fl)

# add phenodata from the spreadsheet
pData(fs) = cbind(pData(fs), tab[1:ninst,])

# create a panoply model
# take default parameters
#    nRecursions = 12 (gives 2^12 or 4096 bins)
#    perplexity = a value controlling tSNE
parms = colnames(fs)[c(7:9, 11:22)]    # all FL except LIVEDEAD
pan = panoply(fcs = fs, parameters = parms, nclust = 50)

# we now have a model.  Roll back through and map each sample to the model
ninst = nrow(tab)
nclust = max(pan$clustering$clst)

sadistics = matrix(NA, nrow = ninst, ncol = nclust)
for (i in 1:ninst) {
  cat("mapping", i, "of", ninst, "...")
  fn = tight(gated_base, tab$FCS.File[i])
  ff = read.FCS(fn)      # read the full flowFrame
  sadistics[i, ] = panoply_map_sample(ff, pan)$fractions
  cat("done.\n")
}

# name the columns
colnames(sadistics) = tight("cluster_", 1:ncol(sadistics))

# make it a dataframe for statistical modeling
sadistics = data.frame(sadistics)

# tack on ages of subjects
# exploit the fact that the rows of the sadistics matrix are in the same
# order as the rows of tab
sadistics = cbind(age = tab$Age, sadistics)

# create a categorical variable, for young (<= 35 years), old (>= 65 years)
# and middle (everything else).
# Prepending the numbers so that the boxplots will sort from young - middle - old
age_cat = rep('2_middle', length = nrow(tab))
age_cat[which(tab$Age <= 35)] = '1_young'
age_cat[which(tab$Age >= 65)] = '3_old'
age_cat = factor(age_cat)

sadistics = cbind(age_range = age_cat, sadistics)

# calculate p-values using the wilcoxon rank sum test
# with Benjamini-Hochberg adjustment for multiple comparisons
idx_young = which(sadistics$age_range == '1_young')
idx_old = which(sadistics$age_range == '3_old')
pval = padj = vector('numeric')
for (i in 1:nclust) {
  clus = tight("cluster_", i)
  pval[i] = wilcox.test(sadistics[idx_young, clus], sadistics[idx_old, clus], exact = FALSE)$p.value
}
padj = p.adjust(pval, method = "BH")

# make some boxplots
opar = par(mfrow = c(8, 7), mar = c(2, 2, 2, 0))

for (i in 1:nclust) {
  clus = tight("cluster_", i)

  bg_col = ifelse(padj[i] <= 0.05, "red", "black")
  fmla = formula(tight(clus, " ~ age_range"))
  tit = sprintf("%s (%.1e)", clus, padj[i])
  boxplot(fmla, data = sadistics, col = c("lightgreen", "lightblue", "pink"), ylab = '', main = tit, col.main = bg_col)
}
par(opar)




