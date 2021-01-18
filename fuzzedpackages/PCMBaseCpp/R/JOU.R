# Copyright 2018 Venelin Mitov
#
# This file is part of PCMBaseCpp.
#
# PCMBase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PCMBase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PCMBase.  If not, see <http://www.gnu.org/licenses/>.

#'@export
PCMInfoCpp.JOU <- function(X, tree, model, 
                           SE = matrix(0.0, PCMNumTraits(model), PCMTreeNumTips(tree)),
                           metaI = PCMInfo(X, tree, model, SE, verbose, preorder=PCMTreePreorderCpp(tree)), 
                           verbose = FALSE, ...) {
  
  if(is.Transformable(model)) {
    model <- PCMApplyTransformation(model)
  }
  
  metaI$pcListInt <- PCListInt(metaI$pc)
  
  res <- c(metaI, cppObject = PCMBaseCpp__QuadraticPolyJOU$new(X, tree, model, metaI))
  
  res$TraverseTree = res$cppObject$TraverseTree
  res$StateAtNode = res$cppObject$StateAtNode
  
  class(res) <- c("PCMInfoCpp", class(metaI))
  res
}

