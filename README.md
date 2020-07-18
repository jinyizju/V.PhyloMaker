V.PhyloMaker

This R package makes phylogenetic hypotheses for a user-specified list of species, by providing multiple ways of binding the species (tips) to a backbone phylogeny, as described in Jin & Qian (2019).

How to install this package? 
This package can be installed in R either using the install_github function in the ‘devtools’ package, or using the githubinstall function in the ‘githubinstall’ package. For example, the R code for installation of V.PhyloMaker using the ‘devtools’ package is as follows:
devtools::install_github("jinyizju/V.PhyloMaker")

Currently, this package includes two major components: (1) the functions, including ‘phylo.maker’, ‘bind.relative’, ‘build.nodes.1’, ‘build.nodes.2’, ‘at.node’, ‘int.node’ and ‘ext.node’; and (2) the data, including 'GBOTB.extended', ‘tips.info’, ‘nodes.info.1’ and ‘nodes.info.2’.

(1) The functions

The main function ‘phylo.maker’ performs the task of binding tips to a backbone phylogeny, which by default is the maga-phylogeny 'GBOTB.extended' of vascular plants embedded in the package.

The function ‘bind.relative’ performs the task of binding tips to the their species- and genus-level closest relative in the backbone phylogeny, when information is provided. ‘bind.relative’ can work alone or with ‘phylo.maker’, for a list of species, as examplified in the function help page and Jin & Qian (2019).

The functions ‘build.nodes.1’ and ‘build.nodes.2’ extract the genus- and family- level node and age information in a phylogeny in different ways. The information extracted are used by ‘phylo.maker’ and ‘bind.relative’ in binding tips to backbone phylogeny.

The functions ‘at.node’, ‘int.node’ and ‘ext.node’ bind a tip to a phylogeny, at different places. These three functions are fast in binding tip to large phylogenies.

(2) The data

The four data sets work with the functions to make phylogenetic hypotheses for a user-specified list of vascular plant species. 

‘GBOTB.extended’ is a mega-tree derived from two recently published mega-trees, described in the help page, and includes 74,531 species and all families of extant vascular plants, is the largest dated phylogeny for vascular plants.

‘tips.info’ is a data frame that contains the information of all the tips, as well as as their genus and family assignments, in ‘GBOTB.extended’.

‘nodes.info.1’ and ‘nodes.info.2’ are the genus- and family-level node and age information in ‘GBOTB.extended’ extracted by ‘build.nodes.1’ and ‘build.nodes.2’, respectively.



Citation:
Jin, Yi, and Hong Qian. 2019. V.PhyloMaker: an R package that can generate very large phylogenies for vascular plants. Ecography 42: 1353–1359. doi: 10.1111/ecog.04434
Smith, Stephen A. and Joseph W. Brown. 2018. Constructing a broadly inclusive seed plant phylogeny. American Journal of Botany 105(3): 302-314. doi: 10.1002/ajb2.1019
Zanne, Amy E., David C. Tank, William K. Cornwell, Jonathan M. Eastman, Stephen A. Smith, Richard G. FitzJohn, Daniel J. McGlinn, Brian C. O’Meara, Angela T. Moles, Peter B. Reich, Dana L. Royer, Douglas E. Soltis, Peter F. Stevens, Mark Westoby, Ian J. Wright, Lonnie Aarssen, Robert I. Bertin, Andre Calaminus, Rafaël Govaerts, Frank Hemmings, Michelle R. Leishman, Jacek Oleksyn, Pamela S. Soltis, Nathan G. Swenson, Laura Warman, Jeremy M. Beaulieu. 2014. Three keys to the radiation of angiosperms into freezing environments. Nature, 506(7486): 89-92. doi: 10.1038/nature12872
