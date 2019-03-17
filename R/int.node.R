int.node <-
function(phylogeny, location.node, tip.label, node.label = NULL, position = 0.5)
{
  phylo <- reorder(phylogeny)
  if (!is.numeric(location.node))
    {
      location.node <- which(phylo$node.label == location.node) + length(phylo$tip.label)
    }
  a <- location.node - length(phylo$tip.label)
  EL <- branching.times(phylo)[a]
  a0 <- a + length(phylo$tip.label)
  a1 <- which(phylo$edge[, 2] == a0)
  aa <- length(which(phylo$edge[1 : a1, 2] <= length(phylo$tip.label)))
  if (is.null(phylo$node.label))
    {
      if (!is.null(node.label))
        {
          nL <- rep(NA, phylo$Nnode + 1)
          nL[a] <- node.label
        }
      if (is.null(node.label))
        {
          nL <- NULL
        }
    }
  if (!is.null(phylo$node.label))
    {
      if (!is.null(node.label)) 
        {
          nL <- c(phylo$node.label[1 : (a - 1)], node.label, phylo$node.label[a : phylo$Nnode])
        }
      if (is.null(node.label)) 
        {
          nL <- c(phylo$node.label[1 : (a - 1)], NA, phylo$node.label[a : phylo$Nnode])
        }
    }
  eG0 <- matrix(c(phylo$edge[a1, 1] + 1, a0 + 1, a0 + 1, aa + 1, a0 + 1, a0 + 2), nrow = 3, byrow = T)
  eG <- matrix(phylo$edge[(a1 + 1) : dim(phylo$edge)[1], ], ncol = 2)
  eG[, 1] <- eG[, 1] + 1
  s <- which(eG[, 1] > a0)
  eG[, 1][s] <- eG[, 1][s] + 1
  s <- which(eG[, 2] > aa)
  eG[, 2][s] <- eG[, 2][s] + 1
  s <- which(eG[, 2] > (a0 + 1))
  eG[, 2][s] <- eG[, 2][s] + 1
  eG <- rbind(eG0, eG)
  eL <- c(phylo$edge.length[a1] * (1 - position), phylo$edge.length[a1] * position + EL, phylo$edge.length[a1] * position, phylo$edge.length[(a1 + 1) : length(phylo$edge.length)])
  tL <- c(tip.label, phylo$tip.label[(aa + 1) : length(phylo$tip.label)])
  if (a1 > 1) 
    { 
      eGn <- matrix(phylo$edge[1 : (a1 - 1), ], ncol = 2)
      eGn[,1] <- eGn[, 1] + 1
      s <- which(eGn[,1] > (a0 + 1))
      eGn[, 1][s] <- eGn[, 1][s] + 1
      s <- which(eGn[, 2] > aa)
      eGn[, 2][s] <- eGn[, 2][s] + 1
      s <- which(eGn[, 2] > (a0 + 1)) 
      eGn[, 2][s] <- eGn[, 2][s] + 1
      eG <- rbind(eGn, eG) 
      eL <- c(phylo$edge.length[1 : (a1 - 1)], eL)
    }
  if (aa > 0)
    {
      tL <- c(phylo$tip.label[1 : aa], tL)
    }
  phylo$edge <- eG
  phylo$tip.label <- tL
  phylo$edge.length <- eL
  phylo$node.label <- nL
  phylo$Nnode <- phylo$Nnode + 1
  return(phylo)
}
