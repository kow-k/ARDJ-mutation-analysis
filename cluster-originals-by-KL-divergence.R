# KL-divergnce clustering on 36 originals

require(philentropy)

##

resp.all.d <- core.data[c("S.ID", "Origin", "Edit", range.names)]
resp.all.d

# select class -----

select.label <- "o"
selected.d <- subset(resp.all.d, Edit == select.label)
selected.d
#
dm <- selected.d[range.names]
# rescale
dm <- dm/rowSums(dm)
#
rownames(dm) <- selected.d$S.ID
dm

## heatmap of all
D <- c()
for (i in 1:nrow(dm)) {
   u <- dm[i, ]
   for (j in 1:nrow(dm)) {
      v <- dm[j, ]
      w <- rbind(u, v)
      kld <- KL(as.matrix(w)) # Crucially
      D <- c(D, kld)
   }
}
if (debugged) print(D)

# convert D to matrix M ----
M <- matrix(D, nrow(dm), nrow(dm))

# give names to M ----
table(selected.d$Edit)
label.text <- paste(selected.d$Edit, selected.d$S.ID, sep = ": ")
label.text

rownames(M) <- label.text
colnames(M) <- label.text
if (verbose) print(M)
#
heatmap(M, cex.main = cex.val, cex.sub = 1.0, cex.lab = 0.7,
        main = sprintf("KL-divergence heatmap of %s stimuli", select.label))


### end of script