
#' Basic summary functions for lvec objects
#'
#' These functions should behave as their regular counterparts. 
#'
#' @param x an \code{\link{lvec}} object
#' @param na.rm logical indicating whether missing values should be ignored
#' @param ... ignored.
#'
#' @rdname ops_summary
#' @import lvec
#' @export
all.lvec <- function(x, ..., na.rm = FALSE) {
  chunkwise(x, 
    init = function(x) {
      TRUE
    }, 
    update = function(state, x) {
      state & all(x, na.rm = na.rm)
    }, 
    final = function(state) {
      state
    }
  )
}

#' @rdname ops_summary
#' @export
any.lvec <- function(x, ..., na.rm = FALSE) {
  chunkwise(x, 
    init = function(x) {
      FALSE
    }, 
    update = function(state, x) {
      state | any(x, na.rm = na.rm)
    }, 
    final = function(state) {
      state
    }
  )
}

#' @rdname ops_summary
#' @export
prod.lvec <- function(x, ..., na.rm = FALSE) {
  chunkwise(x, 
    init = function(x) {
      1.0
    }, 
    update = function(state, x) {
      state * prod(x, na.rm = na.rm)
    }, 
    final = function(state) {
      state
    }
  )
}


#' @rdname ops_summary
#' @export
sum.lvec <- function(x, ..., na.rm = FALSE) {
  chunkwise(x, 
    init = function(x) {
      0
    }, 
    update = function(state, x) {
      state + sum(x, na.rm = na.rm)
    }, 
    final = function(state) {
      state
    }
  )
}

#' @rdname ops_summary
#' @export
mean.lvec <- function(x, ..., na.rm = FALSE) {
  chunkwise(x, 
    init = function(x) {
      c(0.0, 0.0)
    }, 
    update = function(state, x) {
      c(state[1] + sum(x, na.rm = na.rm),
        state[2] + sum(!is.na(x)))
    }, 
    final = function(state) {
      state[1]/state[2]
    }
  )
}

#' @rdname ops_summary
#' @export
max.lvec <- function(x, ..., na.rm = FALSE) {
  chunkwise(x, 
    init = function(x) {
      -Inf
    }, 
    update = function(state, x) {
      max(state, max(x, na.rm = na.rm))
    }, 
    final = function(state) {
      state
    }
  )
}

#' @rdname ops_summary
#' @export
min.lvec <- function(x, ..., na.rm = FALSE) {
  chunkwise(x, 
    init = function(x) {
      Inf
    }, 
    update = function(state, x) {
      min(state, min(x, na.rm = na.rm))
    }, 
    final = function(state) {
      state
    }
  )
}

#' @rdname ops_summary
#' @export
range.lvec <- function(x, ..., na.rm = FALSE) {
  chunkwise(x, 
    init = function(x) {
      c(Inf, -Inf)
    }, 
    update = function(state, x) {
      c(min(state[1], min(x, na.rm = na.rm)),
        max(state[2], max(x, na.rm = na.rm)))
    }, 
    final = function(state) {
      state
    }
  )
}

