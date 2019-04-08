## Replication material for Wagner & Zeileis (2019, German Economic Review):
## Heterogeneity and Spatial Dependence of Regional Growth in the EU:
## A Recursive Partitioning Approach


########################
## Least squares tree ##
########################

## load growth data
data("GrowthNUTS2", package = "lagsarlmtree")

## fit recursive partition of simple base model
library("partykit")
m_tree <- lmtree(ggdpcap ~ gdpcap0 + shgfcf + shsh + shsm |
    gdpcap0 + accessrail + accessroad + capital + regboarder + regcoast + regobj1 + cee + piigs,
  data = GrowthNUTS2, minsize = 12, alpha = 0.05)

## visualization
plot(m_tree, tp_args = list(which = 1, linecol = "red"))

## estimated coefficients
coef(m_tree)

## full summary of estimated models in all terminal nodes
summary(m_tree)

## all parameter stability tests conducted
library("strucchange")
sctest(m_tree)

## construct table version of summary in terminal nodes
cf <- do.call("rbind", lapply(summary(m_tree), function(x) t(x$coefficients[, 1:2])))
rp <- data.frame(
  node = c(3, 5, 6, 7, ""),
  capital = c("no", "no", "no", "yes", ""),
  piigs = c("no", "yes", "yes", "-", ""),
  regcoast = c("-", "no", "yes", "-", ""))[as.vector(rbind(1:4, 5)),]
print(cbind(rp, cf), digits = 2, row.names = FALSE)


######################
## Spatial lag tree ##
######################

## load spatial weights
data("WeightsNUTS2", package = "lagsarlmtree")

## spatial lag model tree
library("lagsarlmtree")
m_stree <- lagsarlmtree(ggdpcap ~ gdpcap0 + shgfcf + shsh + shsm |
    gdpcap0 + accessrail + accessroad + capital + regboarder + regcoast + regobj1 + cee + piigs,
  data = GrowthNUTS2, listw = WeightsNUTS2$invw,
  minsize = 12, alpha = 0.05)
print(m_stree)
plot(m_stree, tp_args = list(which = 1, linecol = "red"))

## all parameter stability tests conducted
sctest(m_stree)

## coefficients after spatial filtering
cf2 <- cbind(coef(m_stree$lagsarlm)[-1], sqrt(diag(vcov(m_stree$lagsarlm)))[-1])
rownames(cf2)[1:4] <- paste0(rownames(cf2)[1:4], ":(Intercept)")
cf2 <- matrix(as.vector(t(cf2)), nrow = 8,
  dimnames = list(
    rep(unique(substr(rownames(cf2), 6, 6)), each = 2),
    unique(substr(rownames(cf2), 8, nchar(rownames(cf2))))
  )
)
rp2 <- data.frame(
  node = c(3, 5, 6, 7, ""),
  capital = c("no", "no", "no", "yes", ""),
  piigs = c("no", "yes", "yes", "-", ""),
  regcoast = c("-", "no", "yes", "-", ""))[as.vector(rbind(1:4, 5)),]
print(cbind(rp2, cf2), digits = 2, row.names = FALSE)


###########################################
## Robustness checks wrt spatial weights ##
###########################################

## set up listw objects for all trees and models
## (exclude listw objects for which observations without neighbors occur)
trees <- models <- Filter(x = WeightsNUTS2,
  f = function(w) all(rowSums(listw2mat(w)) > 0))

## fit all trees
## (tweak alpha so that all splits of interest are carried out)
for(i in names(trees)) trees[[i]] <- lagsarlmtree(
  ggdpcap ~ gdpcap0 + shgfcf + shsh + shsm |
    gdpcap0 + accessrail + accessroad + capital + regboarder + regcoast + regobj1 + cee + piigs,
  data = GrowthNUTS2, listw = trees[[i]], minsize = 12, maxdepth = 4,
  alpha = switch(substr(i, 1, 3), "knn" = 0.10, "inv" = 0.085, 0.07))

## extract p-value information
splitp <- function(lagsarlmtree) {
  n <- lagsarlmtree$tree$node
  p <- c(n$info$p.value, n$kids[[1]]$info$p.value,
    n$kids[[1]]$kids[[1]]$info$test["p.value", "gdpcap0"],
    n$kids[[1]]$kids[[2]]$info$p.value)
  names(p) <- names(lagsarlmtree$tree$data)[c(n$split$varid,
    n$kids[[1]]$split$varid, 2, n$kids[[1]]$kids[[2]]$split$varid)]
  p
}
treep <- t(sapply(trees, splitp))
summary(treep)

## visualize p-values in node 2 (piigs), node 4 (regcoast), and node 3 (gdpcap0)
pplot <- function(x, i, ...) {
  main <- sprintf("Node %s: %s", i, colnames(x)[i])
  x <- x[,i]
  ylim <- c(0, max(0.1, max(x)))
  plot(x, type = "h", ylim = ylim, axes = FALSE, xlab = "", ylab = "p-value", yaxs = "i", main = main, ...)
  text(seq_along(x), par("usr")[3], labels = names(x), srt = 45, adj = c(1.15, 1.15), xpd = TRUE, cex = 0.9)
  axis(2)
  axis(1, at = c(-1000, 1000))
  abline(h = c(0.05, 0.1), lty = 2)
}
par(mfrow = c(3, 1))
pplot(treep, 2)
pplot(treep, 4)
pplot(treep, 3)

## visualize coefficient estimates
## (completely fixing the tree structure)
GrowthNUTS2$.tree <- factor(predict(m_stree, type = "node"))
models <- Filter(x = WeightsNUTS2,
  f = function(w) all(rowSums(listw2mat(w)) > 0))
models <- sapply(models, function(w) coef(lagsarlm(ggdpcap ~ 0 + .tree / (gdpcap0 + shgfcf + shsh + shsm), data = GrowthNUTS2, listw = w)
))
rho <- models[1,]
models <- t(models[-1, ])
colnames(models)[1:4] <- paste0(colnames(models)[1:4], ":(Intercept)")
colnames(models) <- gsub(".tree", "", colnames(models), fixed = TRUE)
ci <- confint(trees$invw$lagsarlm)[-1, ]

k <- ncol(models)/4
par(mfrow = c(1, k), oma = c(0, 5, 0, 0), mar = c(2, 2, 4, 1))
for(i in 1:k) {
  ki <- (i - 1) * 4 + 1:4
  boxplot(models[, ki], xlab = "",
    main = substr(colnames(models)[ki[1]], 3, nchar(colnames(models)[ki[1]])),
    names = substr(colnames(models)[ki], 1, 1), ylim = range(ci[ki,]),
    border = "gray", lwd = 1.5)
  points(1:4, models["invw", ki], pch = 1, cex = 1.3)
  points(1:4, coef(m_tree)[, i], pch = 2, cex = 1.3)
  arrows(1:4, ci[ki, 1], 1:4, ci[ki, 2], code = 3, length = 0.1, angle = 90, col = 1)
  abline(h = 0, col = "black", lty = 2)
  if(i == 1) {
    mtext("Coefficient estimate", side = 2, line = 3)
    legend("topright", c("all W", "inv W", "OLS"),
      pch = c(NA, 1, 2), col = c("gray", "black", "black"), lty = c(1, 1, NA), lwd = 2, bty = "n")
  }
}

