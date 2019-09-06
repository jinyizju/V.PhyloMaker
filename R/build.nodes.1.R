build.nodes.1 <-
function(tree, tips)
{
  tips <- data.frame(as.matrix(tips), stringsAsFactors = FALSE)
  dimnames(tips)[[2]] <- tolower(dimnames(tips)[[2]])
  tips$species <- gsub("(^[[:alpha:]])", "\\U\\1", tips$species, perl = TRUE)
  tips$genus <- gsub("(^[[:alpha:]])", "\\U\\1", tips$genus, perl = TRUE)
  tips$family <- gsub("(^[[:alpha:]])", "\\U\\1", tips$family, perl = TRUE)
  tree$node.label <- paste("N", 1 : tree$Nnode, sep="")
  node.dep <- branching.times(tree)
  tips <- tips[match(tree$tip.label, tips$species), ]
  tips$No. <- as.integer(dimnames(tips)[[1]])
  clustF <- cluster(tips$family)
  tips$numF <- clustF
  clustG <- cluster(tips$genus)
  tips$numG <- clustG
  clustF.size <- table(clustF)
  clustF.size <- data.frame(numF = as.integer(names(clustF.size)), sizeF = as.integer(clustF.size))
  clustG.size <- table(clustG)
  clustG.size <- data.frame(numG = as.integer(names(clustG.size)), sizeG = as.integer(clustG.size))
  tips1 <- merge(tips, clustF.size, all.x=T)
  tips1 <- merge(tips1, clustG.size, all.x=T)
  tips1 <- tips1[, c("No.", "species", "genus", "family", "numF", "numG", "sizeF", "sizeG")]
  xF <- tips1[!duplicated(tips1$numF), ]
  xF <- xF[rev(order(xF$sizeF)), ]
  xF <- xF[!duplicated(xF$family), ]
  xx <- merge(tips1, data.frame(numF = unique(xF$numF)))
  xxG <- table(xx[!duplicated(xx$genus), ]$numF)
  xxG <- data.frame(numF = as.integer(names(xxG)), gen.n = as.integer(xxG))
  xxS <- table(xx$numF)
  xxS <- data.frame(numF = as.integer(names(xxS)), sp.n = as.integer(xxS))
  xxGS <- merge(xxG, xxS)
  xF1 <- xF[xF$sizeF == 1, ]
  if (dim(xF1)[1] > 0)
    {
      Fn1 <- data.frame(level = "F", family = xF1$family, genus = "", rn = "", rn.bl = 0, bn = "", bn.bl = 0, taxa = xF1$species)
      for (i in 1 : dim(Fn1)[1])
        {
          x0 <- match(Fn1$taxa[i], tree$tip.label)
          x1 <- match(x0, tree$edge[, 2])
          x2 <- tree$edge[x1, 1] - length(tree$tip.label)
          Fn1$bn[i] <- tree$node.label[x2]
          Fn1$rn[i] <- tree$node.label[x2]
          Fn1$rn.bl[i] <- tree$edge.length[x1]
          Fn1$bn.bl[i] <- tree$edge.length[x1]
        }
    }
  if (dim(xF1)[1] == 0)
    {
      Fn1 <- NULL
    }
  xF2 <- xF[xF$sizeF > 1, ]
  if (dim(xF2)[1] > 0)
    {
      tipsF <- merge(tips1, data.frame(numF = xF2$numF))
      tF <- drop.tip(tree, setdiff(tree$tip.label, tipsF$species))
      tipsF <- tipsF[match(tF$tip.label, tipsF$species), ]
      tF$tip.label <- tipsF$family
      n1 <- which(!duplicated(tF$tip.label))
      n2 <- length(tF$tip.label) + 1 - which(!duplicated(rev(tF$tip.label)))
      nn <- sort(unique(c(n1, n2)))
      tF <- drop.tip(tF, setdiff(1 : length(tF$tip.label), nn))
      Fn2 <- data.frame(level = "F", family = xF2$family, genus = "", rn = "", rn.bl = 0, bn = "", bn.bl = 0, taxa = "")
      x <- 1 : length(tF$tip.label)
      for (i in 1 : dim(Fn2)[1])
        {
          h <- which(tF$tip.label == Fn2$family[i])
          tt <- drop.tip(tF, setdiff(x, h))
          Fn2$bn[i] <- tt$node.label
          Fn2$bn.bl[i] <- tt$edge.length[1]
          n <- which(tree$node.label == Fn2$bn[i]) + length(tree$tip.label)
          if (n > min(tree$edge[,1]))
            {
              n1 <- tree$edge[which(tree$edge[,2] == n), 1] - length(tree$tip.label)
              Fn2$rn[i] <- tree$node.label[n1]
              Fn2$rn.bl[i] <- node.dep[Fn2$rn[i]]
            }
          if (n == min(tree$edge[,1]))
            {
              Fn2$rn[i] <- Fn2$bn[i]
              Fn2$rn.bl[i] <- Fn2$bn.bl[i]
            }
         }
    }
  if (dim(xF2)[1] == 0)
    {
      xF2 <- NULL
    }
  Fn <- rbind(Fn1, Fn2)
  Fn <- merge(Fn, xF[, c("family", "numF")])
  Fn <- merge(Fn, xxGS)
  Fn <- Fn[, c("level", "family", "genus", "rn", "rn.bl", "bn", "bn.bl", "gen.n", "sp.n", "taxa")]
  xG <- tips1[!duplicated(tips1$numG), ]
  xG <- xG[rev(order(xG$sizeG)), ]
  xG <- xG[!duplicated(xG$genus), ]
  xx <- merge(tips1, data.frame(numG = unique(xG$numG)))
  xxG <- table(xx$numG)
  xxG <- data.frame(numG = as.integer(names(xxG)), gen.n = 1)
  xxS <- table(xx$numG)
  xxS <- data.frame(numG = as.integer(names(xxS)), sp.n = as.integer(xxS))
  xxGS <- merge(xxG, xxS)
  xG1 <- xG[xG$sizeG == 1, ]
  if (dim(xG1)[1] > 0)
    {
      Gn1 <- data.frame(level = "G", family = xG1$family, genus = xG1$genus, rn = "", rn.bl = 0, bn = "", bn.bl = 0, taxa = xG1$species)
      for (i in 1 : dim(Gn1)[1])
        {
          x0 <- match(Gn1$taxa[i], tree$tip.label)
          x1 <- match(x0, tree$edge[, 2])
          x2 <- tree$edge[x1, 1] - length(tree$tip.label)
          Gn1$bn[i] <- tree$node.label[x2]
          Gn1$rn[i] <- tree$node.label[x2]
          Gn1$rn.bl[i] <- tree$edge.length[x1]
         Gn1$bn.bl[i] <- tree$edge.length[x1]
       }
    }
  if (dim(xG1)[1] == 0)
    {
      Gn1 <- NULL
    }
  xG2 <- xG[xG$sizeG > 1, ]
  if (dim(xG2)[1] > 0)
    {
      tipsG <- merge(tips1, data.frame(numG = xG2$numG))
      tG <- drop.tip(tree, setdiff(tree$tip.label, tipsG$species))
      tipsG <- tipsG[match(tG$tip.label, tipsG$species),]
      tG$tip.label <- tipsG$genus
      n1 <- which(!duplicated(tG$tip.label))
      n2 <- length(tG$tip.label) + 1 - which(!duplicated(rev(tG$tip.label)))
      nn <- sort(unique(c(n1, n2)))
      tG <- drop.tip(tG, setdiff(1:length(tG$tip.label), nn))
      Gn2 <- data.frame(level = "G", family = xG2$family, genus = xG2$genus, rn = "", rn.bl = 0, bn = "", bn.bl = 0, taxa = "")
      x<-1 : length(tG$tip.label)
      for (i in 1 : dim(Gn2)[1])
        {
          h <- which(tG$tip.label == Gn2$genus[i])
          tt <- drop.tip(tG, setdiff(x, h))
          Gn2$bn[i] <- tt$node.label
          Gn2$bn.bl[i] <- tt$edge.length[1]
          n <- which(tree$node.label == Gn2$bn[i]) + length(tree$tip.label)
          if (n > min(tree$edge[,1]))
            {
              n1 <- tree$edge[which(tree$edge[,2] == n), 1] - length(tree$tip.label)
              Gn2$rn[i] <- tree$node.label[n1]
              Gn2$rn.bl[i] <- node.dep[Gn2$rn[i]]
            }
          if (n == min(tree$edge[,1]))
            {
              Gn2$rn[i] <- Gn2$bn[i]
              Gn2$rn.bl[i] <- Gn2$bn.bl[i]
            }
        }
    }
  if (dim(xG2)[1] == 0)
    {
      Gn2 <- NULL
    }
  Gn <- rbind(Gn1, Gn2)
  Gn <- merge(Gn, xG[, c("genus", "numG")])
  Gn <- merge(Gn, xxGS)
  Gn <- Gn[, c("level", "family", "genus", "rn", "rn.bl", "bn", "bn.bl", "gen.n", "sp.n", "taxa")]
  nodes <- rbind(Fn, Gn)
  return(nodes)
}
