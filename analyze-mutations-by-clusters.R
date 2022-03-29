# Layering originals 2: Cluster-wise comparison ----

# load packages ----

require(clusternor)
require(dbscan)
require(Rtsne)
require(plotrix)
require(philentropy) # KL divergence

# parameters ----

# viewing
par.old <- par
rm.fn <- "Gill Sans"
ja.fn <- "HiraKakuPro-W3"
show.s.text <- F

#mar.val <- c(4,5,4,5) # this hides the subtitle text
mar.val <- c(5,5,5,5)
if (show.s.text) {
   par(family = eval(ja.fn), xpd = T, mar = mar.val, cex = 0.7)
} else {
   par(family = eval(rm.fn), xpd = T, mar = mar.val, cex = 0.9)
}

# data preparation ----

o.dm <- as.matrix(o.resp.all[range.names]) # Crucially, as data matrix
class(o.dm)
colnames(o.dm)
nrow(o.dm)

# add cluster assignment ----

#o.dm2 <- cbind(o.dm, cluster = o.clustering)

o.dm2 <- o.dm.clustered
#
#grouped.oids <- c("o01", "o11")
o.cluster.ids

## Clusterwise mutation analysis ----

par(mfrow = c(2,2))

xlab.text <- "rating range"
ylab.text <- "response potential [density]"
cex.val <- 1.0
axis.width <- 0.6

## * select cluster id ----

cluster.id <- 1

cluster.filter <- cluster.id == o.dm2[ , which(colnames(o.dm2) == "cluster")]
if (verbose) cluster.filter
grouped.sids <- o.dm2[cluster.filter, ]
if (verbose) grouped.sids
grouped.oids <- subset(sid.oid.map, S.ID %in% rownames(grouped.sids))[['Origin']]
grouped.oids

# set size of group/cluster
g.size <- length(grouped.oids)
#
sub.clust.method <- clust.methods[2]
sub.clust.method
#
tri.columnal <- T
if (tri.columnal) {
   par(mfrow = c(2,3))
} else {
   par(mfrow = c(2,2))
}
#
reorder.by.r01 <- T
eps.val <- 0.15; minPts.val <- 1
#
use.tSNE <- F
repel.labels <- T
show.positioning <- F

# * run -----
i <- 0
for (o.id in grouped.oids) {
   i <- i + 1
   print(sprintf("processing %s [id: %s]", o.id, i))
   print(e <- subset(core.data, Origin == o.id))
   # reorder instances ---
   if (reorder.by.r01) {
      e <- e[order(e$r3x, decreasing = F), ]
      e <- e[order(e$r23, decreasing = F), ]
      e <- e[order(e$r12, decreasing = T), ]
      e <- e[order(e$r01, decreasing = T), ]
   }
   # get mutation labels
   d.mutation <- e$Edit
   # extract core data
   print(d <- data.frame(e[c("S.ID", range.names)], row.names = 1))
   n.cases <- nrow(d)
   kmax.val <- floor(n.cases/2) # for Xmeans
   perplexity.val <- floor((n.cases - 1)/3) # for tSNE
   if (n.cases > 1) { # needs filtering out singletons
      temp1 <- "Profiles of %s and its %s mutations\n[%s/%s in %s cluster %s]"
      main.text <- sprintf(temp1, o.id, (n.cases - 1),
                           i, g.size, clust.method, cluster.id)
      # matplot 1 (stimulus-wise) ----
      matplot(t(d), type = "b", pch = "o",
              ylim = c(0,1), ylab = ylab.text, xlab = xlab.text)
      title(main = main.text, cex.main = cex.val)
      legend("bottomright", pch = "o", col = 1:n.cases,
             cex = 0.7*cex.val, title = "stimuli",
             legend = paste(d.mutation, rownames(d), sep = ": "),
             inset = ifelse(tri.columnal, c(-0.40,-0.15), c(-0.45,-0.15)))
      
      # matplot 2 (range-proportion) ----
      matplot(d, type = "b", pch = "o", ylim = c(0,1), ylab = ylab.text,
              xlab = "stimulus index [after reordering]")
      title(main = sprintf("Proportions of range values by instance\n[%s/%s in %s cluster %s]",
                           i, g.size, clust.method, cluster.id),
            cex.main = cex.val)
      legend("bottomright", pch = "o", col = 1:length(range.names),
             cex = 0.7*cex.val, title = "stimuli",
             legend = range.names,
             inset = ifelse(tri.columnal, c(-0.3,-0.15), c(-0.45,-0.15)))
      
      ## sub clustering ----
      k <- nrow(d)
      dm <- as.matrix(d)
      if (sub.clust.method == "DBSCAN") {
         require(dbscan)
         d.clustered <- dbscan(dm, eps = eps.val, minPts = minPts.val)
         sub.text <- sprintf("clustered by DBSCAN (eps = %s, minPts = %s)",
                 eps.val, minPts.val)
      } else if (sub.clust.method == "Xmeans") {
         require(clusternor)
         d.clustered <- Xmeans(dm, kmax = kmax.val)
         sub.text <- sprintf("clustered by Xmeans (kmax = %s)", kmax.val)
      } else if (sub.clust.method == "FuzzyCMeans") {
         require(clusternor)
         d.clustered <- FuzzyCMeans(dm, centers = k)
         sub.text <- sprintf("clustered by FuzzyCMeans (centers = %s)", kmax.val)
      } else {
         stop("Cluster method unspecified")
      }
      d.cluster.map <- d.clustered$cluster
      if (0 %in% d.cluster.map) {
         d.cluster.map <- d.cluster.map + 1 # handle noises
      }
      print(d.cluster.map)
      n.clust <- length(unique(d.cluster.map))
      ###
      # rotation
      if (use.tSNE) {
         require(Rtsne)
         d.tSNE <- Rtsne(d, check_duplicates = F, pca_scale = T,
                         dims = dims.val, perplexity = perplexity.val)
         d.rot <- d.tSNE$Y
         rownames(d.rot) <- rownames(d)
         colnames(d.rot) <- sprintf("Dim %s", 1:dims.val)
      } else {
         d.pca <- prcomp(d)
         d.rot <- d.pca$x
      }
      if (verbose) print(d.rot)
      ###
      # plot 1: mutation ----
      if (show.positioning) {
         if (use.tSNE) {
            plot(d.rot, col = d.mutation, xpd = T)
         } else {
            plot(d.rot, col = d.mutation, xpd = T,
                 xlim = c(-axis.width, axis.width),
                 ylim = c(-axis.width, axis.width))
         }
         # add title
         temp2 <- "%s of %s and its %s mutations\n[%s/%s in %s cluster %s]"
         main.text <- sprintf(temp2, ifelse(use.tSNE, sprintf("tSNE (pplx: %s)", perplexity.val),
                                     "PCA"),
                              o.id, (n.cases - 1), i, g.size, clust.method, cluster.id)
         title(main = main.text, sub = sub.text, cex.main = cex.val)
         # add legend
         legend("bottomright", title = "mutations",
                cex = 0.7*cex.val, pch = "o", col = d.mutation,
                legend = unique(d.mutation),
                inset = ifelse(tri.columnal, c(-0.3,-0.15), c(-0.35,-0.15))
         )
         # add labels
         label.text <- paste(d.mutation, rownames(d.rot), sep = ": ")
         if (repel.labels) {
            spread.labels(d.rot[ ,1], d.rot[ ,2],
                          linecol = 2, cex = 0.7*cex.val, labels = label.text)
         } else {
            text(d.rot, pos = 3, cex = 0.8*cex.val, labels = label.text)
         }
         }
      ###
      # plot 2: clustering ----
      if (use.tSNE) {
         plot(d.rot, col = d.cluster.map, xpd = T)
      } else {
         plot(d.rot, col = d.cluster.map, xpd = T,
              xlim = c(-axis.width, axis.width),
              ylim = c(-axis.width, axis.width))
      }
      # add title
      temp3 <- "%s of %s and its %s mutations\n[%s/%s in %s cluster %s]"
      main.text <- sprintf(temp3,
                           ifelse(use.tSNE,
                                  sprintf("tSNE (pplx: %s)", perplexity.val),
                                  "PCA"),
                           o.id, (n.cases - 1), i, g.size, clust.method, cluster.id)
      title(main = main.text, sub = sub.text, cex.main = cex.val)
      # add labels
      label.text <- paste(d.mutation, rownames(d.rot), sep = ": ")
      if (repel.labels) {
         spread.labels(d.rot[ ,1], d.rot[ ,2],
                       linecol = 2, cex = 0.7*cex.val, labels = label.text)
      } else {
         text(d.rot, pos = 3, cex = 0.8*cex.val, labels = label.text)
      }
      # add legend
      legend("bottomright", cex = 0.7*cex.val, pch = "o", col = 1:n.clust,
             title = "clusters", legend = sort(unique(d.cluster.map)),
             inset = ifelse(tri.columnal, c(-0.30,-0.15), c(-0.35,-0.15))
      )
      ###
      # violin plot
      if (show.violin.plots) {
         temp5 <- "Violin plot of %s and its variations (Cluster %s)"
         violin_plot(d, col = "lightblue", box_width = 0.1,
               box_col = "blue", median_col = "magenta", xlab = xlab.text,
               main = sprintf(temp5, o.id, cluster.id))
      }
   } else {
      print(sprintf("ignored %s [id: %s] due to lack of cases", o.id, i))
   }
}



## end of script