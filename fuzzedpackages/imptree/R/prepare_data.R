# Copyright (C) 2012 - 2018  Paul Fink
#
# This file is part of imptree.
#
# imptree is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# imptree is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with imptree.  If not, see <https://www.gnu.org/licenses/>.

# Prepare data by given formula or imptree object
# In case of supplied imptree object its formula is used for creating the
# metadata of the Dataset object
prepare_data <- function(object, data, weights, subset, ...) {
  
  # constructing a data.frame according to the supplied formula and na.action
  if(missing(object) || (!inherits(object, "imptree")  && !inherits(object, "formula"))) {
    stop("argument 'object' must be a formula or of class \"imptree\"",
         domain = "R-imptree")
  }
  if(!inherits(data, "data.frame")) {
    stop("argument 'data' must be of class \"data.frame\"", 
         domain ="R-imptree")
  }
  isTree <- inherits(object, "imptree")
  Call <- match.call()
  m <- match(c("data", "subset"), names(Call), 0L)
  mfCall <- Call[c(1L, m)]
  mfCall$formula <- if(isTree) {
    object$formula
  } else {
    object
  }
  mfCall$na.action <- na.fail
  mfCall[[1L]] <- as.name("model.frame")
  mf <- eval(mfCall, parent.frame())
  if (any(attr(attr(mf, "terms"), "order") > 1L)) {
    stop("trees cannot handle interaction terms",
         domain = "R-imptree")
  }
  
  if(!all(sapply(mf, is.factor))) {
    stop("all variables in the resulting model frame must to be of class \"factor\"",
         domain = "R-imptree")
  }
  
  if(missing(weights) || !(length(weights)>0)) {
    wt <- rep(1, nrow(mf))
  } else if(any(weights < 0)) {
    stop("negative weights not allowed", domain = "R-imptree")
  } else {
    wt <- weights
  }

  cpp_data <- as.matrix(data.frame(lapply(mf, as.integer))) - 1
  storage.mode(cpp_data) <- "integer"
  attr(cpp_data, "nlevels") <- sapply(mf, nlevels)
  attr(cpp_data, "labels") <- lapply(mf, levels)
  attr(cpp_data, "classidx") <- 0L
  attr(cpp_data, "wt") <- wt
  
  cpp_data
}
