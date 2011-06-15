# Random Walk Algorithm
#
# Takes as input the matrix of interactions (adjacency matrix) - binary or <=1
# and list of query genes
#
# Produces as output a vector of scores, representing the probability of a
# random walk landing on each genes for a random walk seeded with the query Genes

randomWalk <- function(intM, queryGenes) {
    UseMethod("randomWalk", intM)
}

randomWalk.matrix <- function(intM, queryGenes) {
    if(sum(!queryGenes %in% row.names(intM))>0) {
	stop("queryGenes contains genes not found in intMat")
    }
    Ng <- dim(intM)[1]
    gamma <- 0.7
    # to get the transition matrix for gene network and normalize
    for (i in 1:Ng) {
      intM[,i] <- intM[,i]/sum(intM[,i])
    }
    p0 <- numeric(length=Ng)
    names(p0) <- row.names(intM)

    p0[queryGenes] <- 1
    p0 <- p0/sum(p0)
    res <- .rwr(t(intM),t(p0),gamma)
    return(drop(res))
}


randomWalk.graphNEL <- function(intM, queryGenes) {
    intM <- as(intM, "matrix") # convert to a matrix
    if(sum(!queryGenes %in% row.names(intM))>0) {
	stop("queryGenes contains genes not found in intMat")
    }
    Ng <- dim(intM)[1]
    gamma <- 0.7
    # to get the transition matrix for gene network and normalize
    for (i in 1:Ng) {
      intM[,i] <- intM[,i]/sum(intM[,i])
    }
    p0 <- numeric(length=Ng)
    names(p0) <- row.names(intM)

    p0[queryGenes] <- 1
    p0 <- p0/sum(p0)
    res <- .rwr(t(intM),t(p0),gamma)
    return(drop(res))
}


