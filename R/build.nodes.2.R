build.nodes.2 <-
function(tree, tips)
{
  tips <- data.frame(as.matrix(tips), stringsAsFactors = FALSE)
  tips$species <- gsub("(^[[:alpha:]])", "\\U\\1", tips$species, perl = TRUE)
  tips$genus <- gsub("(^[[:alpha:]])", "\\U\\1", tips$genus, perl = TRUE)
  tips$family <- gsub("(^[[:alpha:]])", "\\U\\1", tips$family, perl = TRUE)
  dimnames(tips)[[2]] <- tolower(dimnames(tips)[[2]])
  tips <- tips[match(tree$tip.label, tips$species), ]
  tree$node.label <- paste("N", 1 : tree$Nnode, sep = "")
  node.dep <- branching.times(tree)
  Fsn <- table(tips$family)
  Fsn <- data.frame(family = names(Fsn), sp.n = as.integer(Fsn))
  Fsn <- merge(Fsn, tips[!duplicated(tips$family), ])
  Fgn <- table(tips[!duplicated(tips$genus), ]$family)
  Fgn <- data.frame(family = names(Fgn), gen.n = as.integer(Fgn))
  Fn <- merge(Fgn, Fsn)
  Fn1 <- Fn[Fn$sp.n == 1, ]
  if (dim(Fn1)[1] > 0)
    {
      Fn1 <- data.frame(level = "F", family = Fn1$family, genus = "", rn = "", rn.bl = 0, bn = "", bn.bl = 0, gen.n = Fn1$gen.n, sp.n = Fn1$sp.n, taxa = Fn1$species, stringsAsFactors = FALSE)
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
  if (dim(Fn1)[1] == 0)
    {
      Fn1 <- NULL
    }
  tF <- tree
  tF$tip.label <- tips$family
  n1 <- which(!duplicated(tF$tip.label))
  n2 <- length(tF$tip.label) + 1 - which(!duplicated(rev(tF$tip.label)))
  nn <- sort(unique(c(n1, n2)))
  tF <- drop.tip(tF, setdiff(1 : length(tF$tip.label), nn))
  Fn2 <- Fn[Fn$sp.n > 1, ]
  if (dim(Fn2)[1] > 0)
    {
      Fn2 <- data.frame(level = "F", family = Fn2$family, genus = "", rn = "", rn.bl = 0, bn = "", bn.bl = 0, gen.n = Fn2$gen.n, sp.n = Fn2$sp.n, taxa = "", stringsAsFactors = FALSE)
      x <- 1 : length(tF$tip.label)
      for (i in 1:dim(Fn2)[1])
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
  if (dim(Fn2)[1] == 0)
    {
      Fn2 <- NULL
    }
  Fn <- rbind(Fn1, Fn2)
  Fn <- Fn[,c("level", "family", "genus", "rn", "rn.bl", "bn", "bn.bl", "gen.n", "sp.n", "taxa")]
  Gsn <- table(tips$genus)
  Gsn <- data.frame(genus = names(Gsn), sp.n = as.integer(Gsn))
  Gsn <- merge(Gsn, tips[!duplicated(tips$genus), ])
  Gsn$gen.n <- 1
  Gn1 <- Gsn[Gsn$sp.n == 1, ]
  if (dim(Gn1)[1] > 0)
    {
      Gn1 <- data.frame(level = "G", family = Gn1$family, genus = Gn1$genus, rn = "", rn.bl = 0, bn = "", bn.bl = 0, gen.n = Gn1$gen.n, sp.n = Gn1$sp.n, taxa = Gn1$species, stringsAsFactors = FALSE)
      for (i in 1 : dim(Gn1)[1])
        {
          x0 <- match(Gn1$taxa[i], tree$tip.label)
          x1 <- match(x0, tree$edge[,2])
          x2 <- tree$edge[x1, 1] - length(tree$tip.label)
          Gn1$bn[i] <- tree$node.label[x2]
          Gn1$rn[i] <- tree$node.label[x2]
          Gn1$rn.bl[i] <- tree$edge.length[x1]
          Gn1$bn.bl[i] <- tree$edge.length[x1]
        }
    }
  if (dim(Gn1)[1] == 0)
    {
      Gn1 <- NULL
    }
  tG <- tree
  tG$tip.label <- tips$genus
  n1 <- which(!duplicated(tG$tip.label))
  n2 <- length(tG$tip.label) + 1 - which(!duplicated(rev(tG$tip.label)))
  nn <- sort(unique(c(n1, n2)))
  tG <- drop.tip(tG, setdiff(1 : length(tG$tip.label), nn))
  Gn2 <- Gsn[Gsn$sp.n > 1, ]
  if (dim(Gn2)[1] > 0)
    {
      Gn2 <- data.frame(level = "G", family = Gn2$family, genus = Gn2$genus, rn = "", rn.bl = 0, bn = "", bn.bl = 0, gen.n = Gn2$gen.n, sp.n = Gn2$sp.n, taxa = "", stringsAsFactors = FALSE)
      x <- 1 : length(tG$tip.label)
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
  if (dim(Gn2)[1] == 0)
    {
      Gn2 <- NULL
    }
  Gn <- rbind(Gn1, Gn2)
  Gn <- Gn[, c("level", "family", "genus", "rn", "rn.bl", "bn", "bn.bl", "gen.n", "sp.n", "taxa")]
  nodes <- rbind(Fn, Gn)
  return(nodes)
}
