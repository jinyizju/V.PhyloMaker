bind.relative <-
function (sp.list, tree = GBOTB.extended, nodes = nodes.info.1,
    output.sp.list = TRUE)
{
    nodesN <- nodes
    treeX <- tree
    # options(scipen = 999)
    sp.list[sapply(sp.list, is.factor)] <- lapply(sp.list[sapply(sp.list,
        is.factor)], as.character)
    if (any(duplicated(sp.list$species))) {
        warning("Duplicated species detected and removed.")
        print(sp.list$species[duplicated(sp.list$species)])
    }
    sp.list <- sp.list[!duplicated(sp.list$species), ]
    sp.list.original <- sp.list
    sp.list$species <- gsub(" ", "_", sp.list$species)
    sp.list$species <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$species,
        perl = TRUE)
    sp.list$genus <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$genus,
        perl = TRUE)
    sp.list$family <- gsub("(^[[:alpha:]])", "\\U\\1", sp.list$family,
        perl = TRUE)
    sp.list$species.relative <- gsub(" ", "_", sp.list$species.relative)
    sp.list$species.relative <- gsub("(^[[:alpha:]])", "\\U\\1",
        sp.list$species.relative, perl = TRUE)
    sp.list$genus.relative <- gsub("(^[[:alpha:]])", "\\U\\1",
        sp.list$genus.relative, perl = TRUE)
    oriN <- tree$node.label
    tree$node.label <- paste("N", 1:length(tree$node.label),
        sep = "")
    add <- sp.list[which(is.na(match(sp.list$species, tree$tip.label))),
        ]
    if (dim(add)[1] == 0 & length(na.omit(match(sp.list$species,
        tree$tip.label))) == 0)
        stop("Incorrect format of species list.")
    if (length(setdiff(sp.list$species, treeX$tip.label)) == 0 & length(na.omit(match(sp.list$species, tree$tip.label))) > 0)
     {
        print("All species in sp.list are present on tree.")
        splis <- sp.list.original
        treeX <- drop.tip(treeX, setdiff(treeX$tip.label, sp.list$species))
        splis$status <- "prune"
        phyloX <- list(phylo = treeX, species.list = splis)
        return(phyloX)
        stop()
     }
    f <- which(sp.list$species %in% tree$tip.label)
    x <- which(sp.list$species.relative %in% tree$tip.label)
    h.sp <- setdiff(x, f)
    fG <- which(sp.list$genus %in% nodes[nodes$level == "G",
        ]$genus)
    xG <- which(sp.list$genus.relative %in% nodes[nodes$level ==
        "G", ]$genus)
    h.gen <- setdiff(xG, fG)
    h.gen <- setdiff(h.gen, union(x, f))
    h <- union(h.sp, h.gen)
    sel.sp <- sp.list[h.sp, ]
    sel.gen <- sp.list[h.gen, ]
    sel.gen <- sel.gen[!duplicated(sel.gen$genus), ]
    if (dim(sel.gen)[1] > 0) {
        for (i in 1:dim(sel.gen)[1]) {
            n <- match(sel.gen$genus.relative[i], nodes$genus)
            x <- length(tree$tip.label) + which(tree$node.label ==
                nodes$bn[n])
            m <- data.frame(level = "G", family = nodes$family[n],
                genus = sel.gen$genus[i], rn = nodes$bn[n], rn.bl = nodes$rn.bl[n],
                bn = nodes$bn[n], bn.bl = nodes$bn.bl[n], gen.n = 1,
                sp.n = 1, taxa = sel.gen$species[i])
            nodesN <- rbind(nodesN, m)
            tree <- at.node(tree, x, sel.gen$species[i])
        }
    }
    if (dim(sel.sp)[1] > 0) {
        for (i in 1:dim(sel.sp)[1]) {
            n <- which(tree$edge[, 2] == match(sel.sp$species.relative[i],
                tree$tip.label))
            tree <- at.node(tree, tree$edge[n, 1], sel.sp$species[i])
        }
    }
    tree$edge.length <- as.numeric(tree$edge.length)
    tree$node.label <- oriN
    status <- rep("", dim(sp.list)[1])
    status[h] <- "add to relative"
    sp.list.original$status.relative <- status
    if (output.sp.list == FALSE)
        sp.list <- NULL
    phylo <- list(phylo = tree, species.list = sp.list.original,
        nodes.info = nodesN)
    phylo[sapply(phylo, is.null)] <- NULL
    return(phylo)
}


