ext.node <-
function(phylogeny, location.tip, tip.label, node.label = NULL, position = 0.5)
{
  phylo <- reorder(phylogeny)
  if (!is.numeric(location.tip))
    {
      location.tip <- which(phylo$tip.label == location.tip)
    }
  a <- location.tip
  a1 <- which(phylo$edge[, 2] == a)
  h <- phylo$edge[a1, 1]         
  if (is.null(phylo$node.label))
    {
      if (!is.null(node.label))
        {
          nL <- rep(NA, phylo$Nnode + 1)
          n <- h - length(phylo$tip.label)
          nL[n + 1] <- node.label
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
          n <- h - length(phylo$tip.label)
          nL <- c(phylo$node.label[1 : n], node.label)
          if (n < phylo$Nnode) 
            {
              nL <- c(nL, phylo$node.label[(n + 1) : phylo$Nnode]) 
            }
        }
      if (is.null(node.label)) 
        {
          n <- h - length(phylo$tip.label)
          nL <- c(phylo$node.label[1 : n], NA)
          if (n < phylo$Nnode) 
            {
              nL <- c(nL,phylo$node.label[(n + 1) : phylo$Nnode]) 
            }
        }
    }
  eG0 <- matrix(c(h + 1, h + 2, h + 2, a, h + 2, a + 1), nrow = 3, byrow = T)
  eG <- matrix(phylo$edge[1 : (a1 - 1), ], ncol = 2)
  eG[, 1] <- eG[, 1] + 1
  s <- which(eG[, 1] > (h + 1))
  eG[, 1][s] <- eG[, 1][s] + 1
  s <- which(eG[, 2] > a)
  eG[, 2][s] <- eG[, 2][s] + 1
  s <- which(eG[, 2] > (h + 1))
  eG[, 2][s] <- eG[, 2][s] + 1
  eG <- rbind(eG, eG0)
  tL <- c(phylo$tip.label[1 : a], tip.label)
  eL <- c(phylo$edge.length[1 : (a1 - 1)],  phylo$edge.length[a1] * (1 - position), phylo$edge.length[a1] * position, phylo$edge.length[a1] * position)
  if (a < length(phylo$tip.label)) 
    { 
      eGn <- matrix(phylo$edge[(a1 + 1) : (dim(phylo$edge)[1]), ], ncol = 2)
      eGn[,1] <- eGn[,1] + 1
      s <- which(eGn[,1] > (h + 1))
      eGn[, 1][s] <- eGn[, 1][s] + 1
      s <- which(eGn[, 2] > a)
      eGn[, 2][s] <- eGn[, 2][s] + 1
      s <- which(eGn[, 2] > (h + 1))
      eGn[, 2][s] <- eGn[, 2][s] + 1
      eG <- rbind(eG, eGn) 
      tL <- c(tL, phylo$tip.label[(a + 1) : length(phylo$tip.label)])
      eL <- c(eL,phylo$edge.length[(a1 + 1) : length(phylo$edge.length)])
    }
  phylo$edge <- eG
  phylo$tip.label <- tL
  phylo$edge.length <- eL
  phylo$Nnode <- phylo$Nnode + 1
  phylo$node.label <- nL
  return(phylo)
}
