  zzphylo.maker <-
  function (sp.list, tree = GBOTB.extended, nodes = nodes.info.1, 
            output.sp.list = TRUE, output.tree = FALSE, scenarios = "S3", 
            r = 1) 
  {
    options(scipen = 999)
    treeX <- tree
    if (is.null(tree$node.labels)) 
      tree$node.label <- rep("", tree$Nnode)
    dimnames(sp.list)[[2]][1:3] <- c("species", "genus", "family")
    sp.list[sapply(sp.list, is.factor)] <- lapply(sp.list[sapply(sp.list, 
                                                                 is.factor)], as.character)
    #return(sp.list)
    sp.list$mergeid = paste( sp.list$family, sp.list$genus, sp.list$species, sep = '_' )
    if (any(duplicated( sp.list$mergeid  ))) {
      print("Duplicated species detected and removed.")
      print(4444)
      print(sp.list$species[duplicated(sp.list$mergeid)])
    }
    #return(0)
    sp.list <- sp.list[!duplicated(sp.list$mergeid), ]
    sp.list.original <- sp.list
    sp.list$species <- gsub(" ", "_", sp.list$species)
    sp.list$species <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$species, 
                            perl = TRUE)
    sp.list$genus <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$genus, 
                          perl = TRUE)
    sp.list$family <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$family, 
                           perl = TRUE)
    rnN <- data.frame(node.label = paste("N", 1:length(tree$node.label), 
                                         sep = ""), oriN = tree$node.label, stringsAsFactors = FALSE)
    nodes[, c("level", "family", "genus", "rn", "bn", "taxa")] <- lapply(nodes[, 
                                                                               c("level", "family", "genus", "rn", "bn", "taxa")], 
                                                                         as.character)
    tree$node.label <- paste("N", 1:length(tree$node.label), 
                             sep = "")
    kk <- c()
    for (i in 1:length(tree$tip.label)) {
      kk <- c(kk, substring(tree$tip.label[i], 1, gregexpr("_", 
                                                           tree$tip.label[i])[[1]][1] - 1))
    }
    m <- data.frame(num = 1:length(kk), genus = kk, species = tree$tip.label)
    m <- merge(m, nodes[, c("genus", "family")])
    mX <- m
    m <- m[, c("genus", "family")]
    m <- m[!duplicated(m$genus), ]
    dimnames(m)[[2]][2] <- "family_in_tree"
    m <- m[, c("genus", "family_in_tree")]
    m0 <- sp.list[!duplicated(sp.list$genus), c("genus", "family")]
    dimnames(m0)[[2]][2] <- "family_in_sp.list"
    mm <- merge(m0, m)
    g <- mm[which(is.na(match(paste(mm$genus, mm$family_in_sp.list, 
                                    sep = "_"), paste(mm$genus, mm$family_in_tree, sep = "_")))), 
            ]
    if (dim(g)[1] > 0) {
      print("Taxonomic classification not consistent between sp.list and tree.")
      print(g)
    }
    add.tip <- sp.list[which(is.na(match(sp.list$species, tree$tip.label))), 
                       ]
    status <- rep("prune", dim(sp.list)[1])
    status[which(is.na(match(sp.list$species, tree$tip.label)))] <- "bind"
    if (dim(add.tip)[1] == 0 & length(na.omit(match(sp.list$species, 
                                                    tree$tip.label))) == 0) 
      stop("Incorrect format of species list.")
    if (length(setdiff(sp.list$species, treeX$tip.label)) == 
        0 & length(na.omit(match(sp.list$species, tree$tip.label))) > 
        0) {
      print("All species in sp.list are present on tree.")
      splis <- sp.list.original
      treeX <- drop.tip(treeX, setdiff(treeX$tip.label, sp.list$mergeid))
      splis$status <- "prune"
      phyloX <- list(scenario.1 = NULL, scenario.2 = NULL, 
                     scenario.3 = NULL, species.list = splis)
      if ("S1" %in% scenarios) {
        phyloX$scenario.1 <- treeX
      }
      if ("S2" %in% scenarios) {
        phyloX$scenario.2 <- treeX
      }
      if ("S3" %in% scenarios) {
        phyloX$scenario.3 <- treeX
      }
      phyloX[sapply(phyloX, is.null)] <- NULL
      return(phyloX)
      stop()
    }
    add.tip$sort <- ""
    add.tip$sort[which(!is.na(match(add.tip$genus, nodes[nodes$level == 
                                                           "G", ]$genus)))] <- "G1"
    add.tip$sort[which(is.na(match(add.tip$genus, nodes[nodes$level == 
                                                          "G", ]$genus)) & !is.na(match(add.tip$family, nodes[nodes$level == 
                                                                                                                "F", ]$family)))] <- "F1"
    add.tip$sort[add.tip$sort == "F1"][duplicated(add.tip[add.tip$sort == 
                                                            "F1", ]$genus)] <- "F2"
    a <- which(add.tip$sort == "")
    if (length(a) > 0) {
      print(paste("Note:", length(a), "taxa fail to be binded to the tree,", 
                  sep = " "))
      print(add.tip$species[a])
      status[match(add.tip$species[a], sp.list$species)] <- "fail to bind"
    }
    sp.list.original$status <- status
    if ("S1" %in% scenarios) {
      t1 <- tree
      rnN1 <- rnN
      nG <- nodes[nodes$level == "G", ]
      nF <- nodes[nodes$level == "F", ]
      data <- add.tip[add.tip$sort == "F1" | add.tip$sort == 
                        "F2", ]
      if (dim(data)[1] > 0) {
        for (i in 1:dim(data)[1]) {
          num <- nF$bn[match(data$family[i], nF$family)]
          t1 <- at.node(t1, location.node = num, tip.label = data$mergeid[i])
        }
      }
      data <- add.tip[add.tip$sort == "G1", ]
      if (dim(data)[1] > 0) {
        for (i in 1:dim(data)[1]) {
          num <- nG$bn[match(data$genus[i], nG$genus)]
          t1 <- at.node(t1, location.node = num, tip.label = data$mergeid[i])
        }
      }
      t1$edge.length <- as.numeric(t1$edge.length)
      tree1 <- t1
      tree1$node.label <- rnN1$oriN[match(tree1$node.label, 
                                          rnN1$node.label)]
      toDrop <- setdiff(1:length(t1$tip.label), which(!is.na(match(t1$tip.label, 
                                                                   sp.list$mergeid))))
      t1 <- drop.tip(t1, tip = toDrop)
      Re <- which(!is.na(match(t1$node.label, rnN1$node.label)))
      noRe <- which(is.na(match(t1$node.label, rnN1$node.label)))
      t1$node.label[Re] <- rnN1$oriN[match(t1$node.label, 
                                           rnN1$node.label)[Re]]
      t1$node.label[noRe] <- ""
    }
    else {
      t1 <- NULL
      tree1 <- NULL
    }
    if ("S2" %in% scenarios) {
      t2r <- replicate(r, list())
      names(t2r) <- paste("run", 1:r, sep = ".")
      tree2r <- replicate(r, list())
      names(tree2r) <- paste("run", 1:r, sep = ".")
      for (o in 1:r) {
        t2 <- tree
        rnN2 <- rnN
        nG <- nodes[nodes$level == "G", ]
        nF <- nodes[nodes$level == "F", ]
        data <- add.tip[add.tip$sort == "F1", ]
        if (dim(data)[1] > 0) {
          for (i in 1:dim(data)[1]) {
            n <- match(data$family[i], nF$family)
            g <- nF$gen.n[n]
            s <- nF$sp.n[n]
            if (g == 1 & s == 1) {
              num <- match(nF$taxa[n], t2$tip.label)
              nlabel <- paste("N", t2$Nnode + 1, sep = "")
              t2 <- ext.node(t2, location.tip = num, tip.label = data$mergeid[i], 
                             node.label = nlabel, position = 2/3)
              nF$gen.n[n] <- g + 1
              nF$sp.n[n] <- s + 1
              x <- which(t2$node.label == nlabel)
              xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel, 
                                                        oriN = ""))
              xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
              rnN2 <- xx
              nF$bn[n] <- nlabel
            }
            else {
              num <- sample(nG$bn[which(nG$family %in% 
                                          data$family[i])], 1)
              t2 <- at.node(t2, location.node = num, tip.label = data$mergeid[i])
              nF$gen.n[n] <- g + 1
              nF$sp.n[n] <- s + 1
            }
          }
        }
        data <- add.tip[add.tip$sort == "F2", ]
        if (dim(data)[1] > 0) {
          for (i in 1:dim(data)[1]) {
            n <- grep(paste(data$genus[i], "_", sep = ""), 
                      t2$tip.label)
            nlabel <- paste("N", t2$Nnode + 1, sep = "")
            if (length(n) == 1) {
              num <- n
              t2 <- ext.node(t2, location.tip = num, tip.label = data$mergeid[i], 
                             node.label = nlabel, position = 1/2)
              x <- which(t2$node.label == nlabel)
              xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel, 
                                                        oriN = ""))
              xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
              rnN2 <- xx
            }
            if (length(n) > 1) {
              num <- t2$edge[which(t2$edge[, 2] == n[1]), 
                             1]
              t2 <- at.node(t2, location.node = num, tip.label = data$mergeid[i])
            }
          }
        }
        data <- add.tip[add.tip$sort == "G1", ]
        if (dim(data)[1] > 0) {
          for (i in 1:dim(data)[1]) {
            n0 <- match(data$genus[i], nG$genus)
            n <- nG$sp.n[n0]
            nlabel <- paste("N", t2$Nnode + 1, sep = "")
            if (n == 1) {
              num <- t2$tip.label[match(nG$taxa[n0], t2$tip.label)]
              t2 <- ext.node(t2, location.tip = num, tip.label = data$mergeid[i], 
                             node.label = nlabel, position = 1/2)
              x <- which(t2$node.label == nlabel)
              xx <- rbind(rnN2[1:(x - 1), ], data.frame(node.label = nlabel, 
                                                        oriN = ""))
              xx <- rbind(xx, rnN2[x:dim(rnN2)[1], ])
              rnN2 <- xx
              nG$sp.n[n0] <- nG$sp.n[n0] + 1
            }
            if (n > 1) {
              num <- which(t2$node.label == nG$bn[n0]) + 
                length(t2$tip.label)
              num1 <- which(t2$edge[, 1] %in% num)
              part1 <- t2$edge[1:min(num1), ]
              n1 <- max(which(part1[, 1] < num), 0) + 
                1
              part2 <- t2$edge[max(num1):dim(t2$edge)[1], 
                               ]
              n2 <- min(which(part2[, 1] < num), dim(part2)[1] + 
                          1) + max(num1) - 2
              sect <- t2$edge[n1:n2, ]
              sect <- sort(unique(c(sect[, 1], sect[, 
                                                    2])))
              sect <- sect[which(sect > length(t2$tip.label))]
              num2 <- sect[sample(1:length(sect), 1)]
              t2 <- at.node(t2, location.node = num2, 
                            tip.label = data$mergeid[i])
            }
          }
        }
        t2$edge.length <- as.numeric(t2$edge.length)
        tree2 <- t2
        tree2$node.label <- rnN2$oriN[match(tree2$node.label, 
                                            rnN2$node.label)]
        toDrop <- setdiff(1:length(t2$tip.label), which(!is.na(match(t2$tip.label, 
                                                                     sp.list$mergeid))))
        t2 <- drop.tip(t2, tip = toDrop)
        Re <- which(!is.na(match(t2$node.label, rnN2$node.label)))
        noRe <- which(is.na(match(t2$node.label, rnN2$node.label)))
        t2$node.label[Re] <- rnN2$oriN[match(t2$node.label, 
                                             rnN2$node.label)[Re]]
        t2$node.label[noRe] <- ""
        t2r[[o]] <- t2
        tree2r[[o]] <- tree2
      }
    }
    else {
      t2r <- NULL
      tree2r <- NULL
    }
    if ("S3" %in% scenarios) {
      t3 <- tree
      rnN3 <- rnN
      nG <- nodes[nodes$level == "G", ]
      nF <- nodes[nodes$level == "F", ]
      data <- add.tip[add.tip$sort == "F1", ]
      if (dim(data)[1] > 0) {
        for (i in 1:dim(data)[1]) {
          n <- match(data$family[i], nF$family)
          g <- nF$gen.n[n]
          s <- nF$sp.n[n]
          if (g == 1 & s == 1) {
            num <- match(nF$taxa[n], t3$tip.label)
            nlabel <- paste("N", t3$Nnode + 1, sep = "")
            t3 <- ext.node(t3, location.tip = num, tip.label = data$mergeid[i], 
                           node.label = nlabel, position = 2/3)
            nF$gen.n[n] <- g + 1
            nF$sp.n[n] <- s + 1
            x <- which(t3$node.label == nlabel)
            xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel, 
                                                      oriN = ""))
            xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
            rnN3 <- xx
            nF$bn[n] <- nlabel
          }
          if (g == 1 & s > 1) {
            nlabel <- paste("N", t3$Nnode + 1, sep = "")
            if ((2/3) * nF$rn.bl[n] <= nF$bn.bl[n]) {
              len <- (nF$rn.bl[n] - nF$bn.bl[n])/2
            }
            if ((2/3) * nF$rn.bl[n] > nF$bn.bl[n]) {
              len <- nF$rn.bl[n] * 2/3 - nF$bn.bl[n]
            }
            port <- len/(nF$rn.bl[n] - nF$bn.bl[n])
            t3 <- int.node(t3, location.node = nF$bn[n], 
                           tip.label = data$mergeid[i], node.label = nlabel, 
                           position = port)
            nF$gen.n[n] <- g + 1
            nF$sp.n[n] <- s + 1
            x <- which(t3$node.label == nlabel)
            xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel, 
                                                      oriN = ""))
            xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
            rnN3 <- xx
            nF$bn[n] <- nlabel
          }
          if (g > 1) {
            t3 <- at.node(t3, location.node = nF$bn[n], 
                          tip.label = data$mergeid[i])
          }
        }
      }
      data <- add.tip[add.tip$sort == "F2", ]
      if (dim(data)[1] > 0) {
        for (i in 1:dim(data)[1]) {
          n <- grep(paste(data$genus[i], "_", sep = ""), 
                    t3$tip.label)
          nlabel <- paste("N", t3$Nnode + 1, sep = "")
          if (length(n) == 1) {
            t3 <- ext.node(t3, location.tip = n, tip.label = data$mergeid[i], 
                           node.label = nlabel, position = 1/2)
            nG$sp.n[match(data$genus[i], nG$genus)] <- length(n) + 
              1
            x <- which(t3$node.label == nlabel)
            xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel, 
                                                      oriN = ""))
            xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
            rnN3 <- xx
            nG$bn[match(data$genus[i], nG$genus)] <- nlabel
          }
          if (length(n) > 1) {
            num <- ancestor(t3, min(n), max(n))
            t3 <- at.node(t3, location.node = num, tip.label = data$mergeid[i])
          }
        }
      }
      data <- add.tip[add.tip$sort == "G1", ]
      if (dim(data)[1] > 0) {
        for (i in 1:dim(data)[1]) {
          n0 <- match(data$genus[i], nG$genus)
          s <- nG$sp.n[n0]
          nlabel <- paste("N", t3$Nnode + 1, sep = "")
          if (s == 1) {
            num <- t3$tip.label[match(nG$taxa[n0], t3$tip.label)]
            t3 <- ext.node(t3, location.tip = num, tip.label = data$mergeid[i], 
                           node.label = nlabel, position = 1/2)
            nG$sp.n[n0] <- nG$sp.n[n0] + 1
            x <- which(t3$node.label == nlabel)
            xx <- rbind(rnN3[1:(x - 1), ], data.frame(node.label = nlabel, 
                                                      oriN = ""))
            xx <- rbind(xx, rnN3[x:dim(rnN3)[1], ])
            rnN3 <- xx
            nG$bn[n0] <- nlabel
          }
          if (s > 1) {
            t3 <- at.node(t3, location.node = nG$bn[n0], 
                          tip.label = data$mergeid[i])
          }
        }
      }
      t3$edge.length <- as.numeric(t3$edge.length)
      tree3 <- t3
      tree3$node.label <- rnN3$oriN[match(tree3$node.label, 
                                          rnN3$node.label)]
      toDrop <- setdiff(1:length(t3$tip.label), which(!is.na(match(t3$tip.label, 
                                                                   sp.list$mergeid))))
      t3 <- drop.tip(t3, tip = toDrop)
      Re <- which(!is.na(match(t3$node.label, rnN3$node.label)))
      noRe <- which(is.na(match(t3$node.label, rnN3$node.label)))
      t3$node.label[Re] <- rnN3$oriN[match(t3$node.label, 
                                           rnN3$node.label)[Re]]
      t3$node.label[noRe] <- ""
    }
    else {
      t3 <- NULL
      tree3 <- t3
    }
    if (output.sp.list == FALSE) {
      sp.list.original <- NULL
    }
    if (output.tree == FALSE) {
      tree1 <- tree2r <- tree3 <- NULL
    }
    phylo <- list(scenario.1 = t1, scenario.2 = t2r, scenario.3 = t3, 
                  species.list = sp.list.original, tree.scenario.1 = tree1, 
                  tree.scenario.2 = tree2r, tree.scenario.3 = tree3)
    phylo[sapply(phylo, is.null)] <- NULL
    return(phylo)
  }
