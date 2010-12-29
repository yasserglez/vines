# vines: R package for multivariate dependence modeling with vines
# Copyright (C) 2010 Yasser González-Fernández <ygonzalezfernandez@gmail.com>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details.
#
# You should have received a copy of the GNU General Public License along with 
# this program. If not, see <http://www.gnu.org/licenses/>.

layoutDVine <- function (vine, graph, vertexTrees, vertexPositions) {
  layoutMatrix <- matrix(nrow = 0, ncol = 2)

  for (j in seq(length = vine@trees)) {
    for (i in seq(length = vine@dimension - j + 1)) {
      layoutMatrix <- rbind(layoutMatrix, c(i - 1, vine@trees - j))
    }
  }

  layoutMatrix
}

plotDVine <- function (x) {
  if (length(x@dimensionNames) > 0) {
    dimNames <- x@dimensionNames
  } else {
    dimNames <- as.character(seq(length = x@dimension))
  }

  graph <- graph.empty(0, directed = FALSE)
  for (j in seq(length = x@trees)) {
    currentCopulaLabels <- character(0)
    for (i in seq(length = x@dimension - j)) {
      conditioned <- paste(dimNames[i], dimNames[i + j], sep = ",")
      conditioning <- paste(dimNames[seq(from = i + 1, to = i + j - 1)], collapse = ",")
      copulaLabel <- paste(conditioned,
          if (j > 1) paste("|", conditioning, sep = "")
              else character(0), 
          sep = "")
      currentCopulaLabels <- c(currentCopulaLabels, copulaLabel)
    }
    
    if (j == 1) {
      vertexLabels <- dimNames
      edgeLabels <- currentCopulaLabels
    } else {
      vertexLabels <- c(vertexLabels, previousCopulaLabels)
      edgeLabels <- c(edgeLabels, currentCopulaLabels)
    }
    previousCopulaLabels <- currentCopulaLabels

    graph <- graph.disjoint.union(graph, 
        graph.tree(x@dimension - j + 1, 1, mode = "undirected"))
  }

  plot(graph,
      vertex.color = "white",
      vertex.label = vertexLabels,
      vertex.label.cex = 1.0,
      vertex.label.color = "black",
      vertex.shape = "rectangle",
      edge.color = "lightgray",
      edge.label = edgeLabels,
      edge.label.cex = 1.0,
      edge.label.color = "black",
      layout = layoutDVine(vine, graph, vertexTrees, vertexPositions))
}


layoutCVine <- function (graph) {
  
}

plotCVine <- function (x) {
  
}
