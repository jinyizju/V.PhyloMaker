at.node <-
function(phylogeny, location.node, tip.label)
{
  phylo <- reorder(phylogeny)
  if (!is.numeric(location.node))
    location.node <- which(phylo$node.label == location.node) + length(phylo$tip.label)
  a <- location.node - length(phylo$tip.label)
  EL <- branching.times(phylo)[a]
  a0 <- a + length(phylo$tip.label)
  a1 <- which(phylo$edge[,1] == a0)[1]
  aa <- length(which(phylo$edge[1 : (a1 - 1), 2] <= length(phylo$tip.label)))
  eG0 <- matrix(c(a0 + 1, aa + 1), nrow = 1)
  eG <- matrix(phylo$edge[a1 : dim(phylo$edge)[1], ], ncol = 2)
  eG[, 1] <- eG[, 1] + 1
  s <- which(eG[, 2] > aa)
  eG[, 2][s] <- eG[, 2][s] + 1
  eG <- rbind(eG0, eG)    
  eL <- c(EL, phylo$edge.length[a1 : length(phylo$edge.length)])
  tL <- c(tip.label, phylo$tip.label[(aa + 1) : length(phylo$tip.label)])
  if (a1 > 1) 
    { 
      eGn <- matrix(phylo$edge[1 : (a1 - 1), ], ncol = 2)
      eGn[, 1] <- eGn[, 1] + 1
      s <- which(eGn[, 2] > aa)
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
  return(phylo)
}
