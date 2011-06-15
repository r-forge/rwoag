\name{randomWalk}
\alias{randomWalk}
\alias{randomWalk.matrix}
\alias{randomWalk.graphNEL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs a Random Walk on a Graph
}
\description{
Performs random walk on a graph, with a restart probability of 0.7
}
\usage{
randomWalk(intM, queryGenes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intM}{
    adjacency matrix of the graph, or object of graph format containing this matrix
}
  \item{queryGenes}{
    genes with which to seed the random walk
}
}
\value{
    Returns a named vector of visitation probabilities for all nodes in the graph
}
\references{
    http://www.ncbi.nlm.nih.gov/pubmed/18371930
}
\author{
    James Perkins
    Ian Morilla

    Maintainer: James Perkins <jperkins@biochem.ucl.ac.uk>
}

\examples{
    # using a graphNEL object
    set.seed(123)
    V <- letters[1:10]
    M <- 1:4
    g1 <- randomGraph(V, M, 0.5)
    visProbs <- randomWalk(g1, V[1:3])

    # using a graph object
    require(graph)
    zz <- rbind(c(0,1,0,1), c(1,0,0,0), c(0,0,0,1), c(1,0,1,0))
    colnames(zz) = rownames(zz) = letters[1:4]
    visProbs <- randomWalk(zz, V[1])
}

\keyword{ graphs }
