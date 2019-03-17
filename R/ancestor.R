ancestor <-
function(tree, tip1, tip2)
{
  ance1 <- tree$edge[which(tree$edge[, 2] == tip1), 1]
  ance2 <- tree$edge[which(tree$edge[, 2] == tip2), 1]
  mi <- min(tree$edge[, 1])
  for (i in 1 : (tree$Nnode - 1))
    {
       if (ance1[i] > mi)  ance1[i + 1] <- tree$edge[which(tree$edge[,2] == ance1[i]), 1]
       else  break()
     }
  for (j in 1 : (tree$Nnode - 1))
    {
       if (ance2[j] > mi)  ance2[j + 1] <- tree$edge[which(tree$edge[,2] == ance2[j]), 1]
       else  break()
    }
  node <- max(intersect(ance1, ance2))
  return(node)
}
