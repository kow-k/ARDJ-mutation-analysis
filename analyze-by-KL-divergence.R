##

debugged <- F

## KL-divergence ----

require(philentropy) # provides KL(..)
#require(transport)  # provides aha(..)

if (cluster.id %in% c(1,2,3)) {
   par(mar = c(5,5,5,5), cex = 0.9)
   # * run ----
   i <- 0
   for (o.id in grouped.oids) {
      i <- i + 1
      print(sprintf("processing %s [id: %s]", o.id, i))
      print(e <- subset(core.data, Origin == o.id))
      # get mutation labels
      d.mutation <- e$Edit
      # extract core data
      print(d <- data.frame(e[c("S.ID", range.names)], row.names = 1))
      n.cases <- nrow(d)
      if (n.cases > 1) { # needs filtering out singletons
         d <- d/rowSums(d) # rescaled for KL-divergence compatibility
         # define labels
         label.text <- paste(d.mutation, rownames(d), sep = ": ")
         
         # define titles
         main.text <- sprintf("KL-divergence cross-table with %s group [%s cluster %s]",
                              o.id, clust.method, cluster.id)
         #
         D <- c()
         OT <- c()
         for (i in 1:nrow(d)) {
            u <- d[i, ]
            for (j in 1:nrow(d)) {
               v <- d[j, ]
               w <- rbind(u, v)
               #apply(w, 2, rowSums)
               kld <- KL(as.matrix(w)) # Crucially
               D <- c(D, kld)
               # OT
               otd <- ada(u, v, C) #
               OT <- c(OT, otd)
            }
         }
         if (debugged) print(D)
         M <- matrix(D, nrow(d), nrow(d))
         #rownames(M) <- rownames(d)
         #colnames(M) <- rownames(d)
         rownames(M) <- label.text
         colnames(M) <- label.text
         print(sprintf("KL-divergence cross-table for %s", o.id))
         print(M)
         heatmap(M, cex.main = cex.val, cex.sub = 1.0, cex.lab = 0.7,
                 main = main.text, sub = sub.text)
      }
   }
   #   
}

### end of script