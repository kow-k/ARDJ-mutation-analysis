## ARDJ mutation analysis
# created on 2021/04/30
# by Kow Kuroda
# modified on 2022/03/29 to make ready for public release

# parameters ----
verbose <- F
debugged <- F
if (debugged) verbose <- T
refresh.data <- T
# select source 
use.s2w <- F # s2u needs to be used, s2w data lacks Origin var

# processing
use.sorted <- T
sampled <- F
sample.n <- 8

# source selection
#use.halves <- T
use.u.half <- T
use.l.half <- F
if (use.l.half) use.u.half <- F
use.all <- F
if (use.all) {
   use.u.half <- F
   use.l.half <- F
} else { use.halves <- T }
#
if (use.halves) {
   use.all <- F
   if (use.l.half) use.u.half <- F
}

# clustering method
clust.methods <- c("Xmeans", "DBSCAN", "FuzzyCMeans")
clust.selector <- 2
clust.method <- clust.methods[clust.selector]
clust.method
#
if (clust.method == "DBSCAN") {
   use.DBSCAN <- T
} else if (clust.method == "Xmeans") {
   use.Xmeans <- T # use with caution, crash-prone
}
clust.method

#
plot.halves <- F
show.s.text <- F
show.violin.plots <- F
repel.labels <- T

# viewing
par.old <- par
rm.fn <- "Gill Sans"
ja.fn <- "HiraKakuPro-W3"
#
par.old <- par
mfrow.val <- c(2,2)
#mar.val <- c(4,5,4,5) # this hides the subtitle text
mar.val <- c(5,5,5,5)
if (show.s.text) {
   par(family = eval(ja.fn), xpd = T, mar = mar.val, cex = 0.7)
} else {
   par(family = eval(rm.fn), xpd = T, mar = mar.val)
}
#
paned <- T
if (paned) {
   par(mfrow = mfrow.val)
   cex.val <- 0.7
}

# graphic parameters -----
xlab.text <- "rating ranges"
ylab.text <- "response potential [density]"


# load packages ----
require(gdata) # provides read.xls
#require(MASS)
require(plotrix) # provides violin_plot and spread.labels
require(dbscan) # provides dbscan(..)
require(Rtsne) # provides tSNE(..)
require(clusternor) # provides Xmeans(..), FuzzyCMeans(..) 
#require(KernSmooth)
#require(FactoMineR)

# Load data ----
#data.dir <- "/Users/kowk/Dropbox/ARDJ/results/responses"

# set target files ----
#if (use.s2w) {
#   file.name <- paste(data.dir, "response-base-s2w.xlsx", sep = "/")
#   sheet.name <- "s2w.r1.filtered"
#} else {
#   file.name <- paste(data.dir, "response-base-s2u.xlsx", sep = "/")
#   sheet.name <- "s2u.sd.filtered1"
#}
#file.name
#sheet.name

file.name <- "./data-s2u-sd-filtered1.csv"
file.name

# reload file ----
if (refresh.data) {
   #resp.data.raw <- read.xls(file.name, sheet = sheet.name, header = T)
   read.csv(file.name)
}
if (verbose) View(resp.data.raw)
nrow(resp.data.raw)

# generate clean data 
resp.data.clean <- subset(resp.data.raw, !(S.ID==""))
if (verbose) View(resp.data.clean)
nrow(resp.data.clean)

# alias ----
resp.data <- resp.data.clean
View(resp.data)

# define variables ----
var.selector <- "S.ID|V.ID|S.TEXT|r01|r12|r23|r3x|Edit|Origin"
range.names <- c("r01", "r12", "r23", "r3x")
mutation.order <- c("o", "s", "n", "v", "p")
#
core.data.init <- resp.data[grepl(var.selector, colnames(resp.data))]
if (verbose) View(core.data.init)

# remove irrelevant rows ----
sprintf("core.data.init contains: %s", 
   nrow(core.data.init)) # returns 300

core.data <- subset(core.data.init, !Origin == "")
sprintf("core.data contains: %s",
   nrow(core.data)) # returns 270


## all included ----

all.dm <- data.frame(core.data[c("S.ID", range.names)], row.names = 1)
str(all.dm)

k.all <- nrow(all.dm)
k.all
all.centers.val <- floor(k.all / 3)
all.centers.val
all.perplexity.val <- floor((k.all - 1)/3)
all.perplexity.val

# * select dbscan paramters -----
eps.val <- 0.07; minPts.val <- 2 # eps = 0.15 turns out be too large, eps = 0.06 too small

# run
if (clust.method == "DBSCAN") {
   all.clustered <- dbscan(all.dm, eps = eps.val, minPts = minPts.val)
   sub.text <- sprintf("clustered by DBSCAN (eps = %s, minPts = %s)",
                       eps.val, minPts.val)
   if (0 %in% unique(all.clustered$cluster)) {
      all.clustered$cluster <- (all.clustered$cluster + 1) # handling noises
   }
} else if (clust.method == "Xmeans") {
   all.clustered <- Xmeans(all.dm, kmax = all.centers.val)
   sub.text <- sprintf("clustered by Xmeans (kmax = %s)", all.centers.val)
} else if (clust.method == "FuzzyCMeans") {
   all.clustered <- FuzzyCMeans(all.dm, centers = all.centers.val)
   sub.text <- sprintf("clustered by FuzzyCMeans (centers = %s)", all.centers.val)
} else { stop("Clustering method unspecified") }
#
all.cluster.map <- all.clustered$cluster
clust.method
all.cluster.map
sort(unique(all.cluster.map))
#
all.cluster.ids <- sort(unique(all.cluster.map))
all.cluster.ids

# assign cluster -----
all.dm.clustered <- cbind(all.dm, cluster = all.cluster.map)


## Analysis using PCA, tSNE ----

# PCA ----

all.pca <- prcomp(all.dm, scale = T)
all.pca
all.rot.pca <- all.pca$x
all.x <- abs(all.pca$rotation)
#x
#colSums(x)
#all.loadings <- sweep(all.x, 2, colSums(all.x), "/")
#all.loadings

# plot PCA
par(mfrow = c(1,1))

plot(all.rot.pca, 
     #xlim = c(-0.5,0.5), ylim = c(-0.5,0.5),
     col = all.cluster.map)
# add title
title(main = sprintf("PCA of all %s stimuli", nrow(all.dm)),
      sub = sub.text)
# add legend
legend("topright", pch = "o", cex = 0.8,
       col = unique(all.cluster.map), title = "clusters",
       legend = sort(unique(all.cluster.map)), inset = c(-0.15,0.05)
)
# add labels
repel.labels <- F
if (repel.labels) {
   spread.labels(all.rot.pca[ ,1], all.rot.pca[ ,2], xpd = T, cex = 0.7,
                 labels = rownames(all.rot.pca))
} else {
   text(all.rot.pca, pos = 3, cex = 0.5, col = "grey",
        labels = rownames(all.rot.pca))
}

# tSNE
all.n.cases <- (nrow(all.dm))
all.n.cases
all.dims.val <- 3
all.perplexity.val <- floor((all.n.cases - 1)/3)
all.perplexity.val
all.tSNE <- Rtsne(all.dm, check_duplicates = F,
                perplexity = all.perplexity.val,
                dim = all.dims.val
                #, pca_scale = T, pca_center = F
)
all.rot.tSNE <- all.tSNE$Y
colnames(all.rot.tSNE) <- sprintf("Dim %s", 1:all.dims.val)
head(all.rot.tSNE)

# plot
plot(all.rot.tSNE, #xlim = c(-0.5,0.5), ylim = c(-0.5,0.5),
     col = all.cluster.map, cex = 1.2)
# add title
title(main = sprintf("tSNE (pplx: %s) of all %s stimuli",
                     all.perplexity.val, all.n.cases),
      sub = sub.text)
# add legend
legend("topright", pch = "o", cex = 0.8,
       col = 1:length(unique(all.cluster.map)),
       title = "clusters", legend = sort(unique(all.cluster.map)),
       inset = c(-0.15,0.05)
)
# add labels
repel.labels.super <- F
if (repel.labels.super) {
   spread.labels(all.rot.tSNE[ ,1], all.rot.tSNE[ ,2], cex = 0.7, xpd = T,
                 labels = rownames(all.dm))
} else {
   text(all.rot.tSNE, pos = 3, cex = 0.6, col = "grey", labels= rownames(all.dm))
}


# sort data by Edit/mutation ----
if (use.sorted) {
   core.data.sorted <- core.data[order(match(core.data$Edit,
                                             mutation.order)), ]
   core.data <- core.data.sorted
}
head(core.data)

# get originals -----
o.only <- subset(core.data, Edit == "o")
nrow(o.only)
head(o.only)
#
o.sid.map <- o.only[c("S.ID", "Origin")]
nrow(o.sid.map)
head(o.sid.map)

# get origin-sid map -----
sid.oid.map <- core.data[c("S.ID", "Origin")]
nrow(sid.oid.map)
head(sid.oid.map)

oids.all <- sort(unique(sid.oid.map$Origin))
oids.all

o.resp.all <- data.frame(o.only[c("S.ID", "Origin", range.names)],
                         row.names = 1)
nrow(o.resp.all)
head(o.resp.all)
if (verbose) View(o.resp.all)

# remove duplicates ----
duplicated.index <- duplicated(o.resp.all)
rownames(o.resp.all)[duplicated.index]
# o.resp.all.clean <- o.resp.all[-duplicated.index, ] # not working
o.resp.all.clean <- o.resp.all[!duplicated.index, ]
nrow(o.resp.all.clean)
duplicated(o.resp.all.clean)

# reorder data by r01 values ----
if (use.sorted) o.resp.all <- o.resp.all[order(o.resp.all$r01,
                                               decreasing = T), ]

# all originals -----
n.origins <- nrow(o.resp.all)
sprintf("number of origins: %s", n.origins)

View(o.resp.all)

# Clustering original ----

## data preparation ----

o.dm <- as.matrix(o.resp.all[range.names]) # Crucially, as data matrix
class(o.dm)
colnames(o.dm)
nrow(o.dm)

## main ----

clust.selector <- 2
clust.method <- clust.methods[clust.selector]
clust.method

k <- nrow(o.dm)
centers.val <- floor(k / 3)
perplexity.val <- floor((k - 1)/3)
# Note: eps = 0.2 turns out to be too large, eps = 0.11 too small
eps.val <- 0.13; minPts.val <- 2
# run
if (clust.method == "DBSCAN") {
   o.clustered <- dbscan(o.dm, eps = eps.val, minPts = minPts.val)
   sub.text <- sprintf("clustered by DBSCAN (eps = %s, minPts = %s)",
                       eps.val, minPts.val)
   if (0 %in% unique(o.clustered$cluster)) {
      o.clustered$cluster <- (o.clustered$cluster + 1) # handling noises
   }
} else if (clust.method == "Xmeans") {
   o.clustered <- Xmeans(o.dm, kmax = k - 1)
   sub.text <- sprintf("clustered by Xmeans (kmax = %s)", k - 1)
} else if (clust.method == "FuzzyCMeans") {
   o.clustered <- FuzzyCMeans(o.dm, centers = k)
   sub.text <- sprintf("clustered by FuzzyCMeans (centers = %s)", centers.val)
} else { stop("Clustering method unspecified") }
#
o.cluster.map <- o.clustered$cluster
clust.method
o.cluster.map
#
o.cluster.ids <- sort(unique(o.cluster.map))
o.cluster.ids

# assign cluster -----
o.dm.clustered <- cbind(o.dm, cluster = o.cluster.map)

## Analysis using PCA, tSNE ----

# PCA ----

o.pca <- prcomp(o.dm, scale = T)
o.pca
o.rot.pca <- o.pca$x
x <- abs(o.pca$rotation)
#x
o.loadings <- sweep(x, 2, colSums(x), "/")
o.loadings

# plot PCA
par(mfrow = c(1,1))

plot(o.rot.pca, 
     #xlim = c(-0.5,0.5), ylim = c(-0.5,0.5),
     col = o.cluster.map)
# add title
title(main = sprintf("PCA of %s originals", k),
      sub = sub.text)
# add legend
legend("topright", pch = "o", cex = 0.8,
       col = unique(o.cluster.map), title = "clusters",
       legend = sort(unique(o.cluster.map)), inset = c(-0.15,0.05)
)
# add labels
repel.labels <- F
if (repel.labels) {
   spread.labels(o.rot.pca[ ,1], o.rot.pca[ ,2], xpd = T,
                 labels = rownames(o.rot.pca))
} else {
   text(o.rot.pca, pos = 3, cex = 0.7, labels = rownames(o.rot.pca))
}

# tSNE ----
n.cases <- (nrow(o.dm))
n.cases
dims.val <- 3
perplexity.val <- floor((n.cases - 1)/3)
perplexity.val
o.tSNE <- Rtsne(o.dm, check_duplicates = F,
                perplexity = perplexity.val,
                dim = dims.val #, pca_scale = T, pca_center = F
)
o.rot.tSNE <- o.tSNE$Y
colnames(o.rot.tSNE) <- sprintf("Dim %s", 1:dims.val)
head(o.rot.tSNE)

# plot
plot(o.rot.tSNE, #xlim = c(-0.5,0.5), ylim = c(-0.5,0.5),
     col = o.cluster.map, cex = 1.2)
# add title
title(main = sprintf("tSNE (pplx: %s) of %s originals", perplexity.val, k),
      sub = sub.text)
# add legend
legend("topright", pch = "o", cex = 0.8,
       col = 1:length(unique(o.cluster.map)),
       title = "clusters", legend = sort(unique(o.cluster.map)),
       inset = c(-0.15,0.05)
)
# add labels
repel.labels.super <- F
if (repel.labels.super) {
   spread.labels(o.rot.tSNE[ ,1], o.rot.tSNE[ ,2],
                 labels = rownames(o.dm), xpd = T)
} else {
   text(o.rot.tSNE, pos = 3, labels= rownames(o.dm), cex = 0.6)
}


## plot profiles clusterwise -----

o.dm.clustered
#o.dm.clustered2 <- cbind(o.dm.clustered, Origin = o.sid.map[, 2])
#head(o.dm.clustered2)

par(mfrow = c(2,2))

cex.val <- 1.0

for (cluster.id in o.cluster.ids) {
   print(sprintf("processing %s cluster: %s", clust.method, cluster.id))
   ds <- subset(as.data.frame(o.dm.clustered), cluster == cluster.id)
   print(ds)
   d <- ds[range.names]
   c.size <- nrow(d)
   # matplot
   main.text <- sprintf("Profiles of %s originals\nin %s cluster %s",
                        c.size, clust.method, cluster.id)
   matplot(t(d), type = "b", pch = "o", cex.main = cex.val, ylim = c(0,1),
           xlab = xlab.text, ylab = ylab.text, main = main.text, sub = sub.text)
   # add legend
   legend.text <- sprintf("%s = %s", rownames(d),
                          o.sid.map[(o.sid.map[ ,'S.ID'] %in% rownames(d)), 'Origin'] )
   legend("topright", title = "stimuli", cex = 0.7, pch = "o", col = 1:c.size,
          legend = legend.text, inset = c(-0.55,-0.2))
}

### end of script
