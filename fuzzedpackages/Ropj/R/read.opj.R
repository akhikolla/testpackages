# We need to preserve the "commands" and "comment" attributes set by the
# C++ code.  Since lapply and as.data.frame don't do that, we wrap them
# and set them again on the returned value.
.preserve.attributes <- function(f) function(x, ...)
	do.call(structure, c(list(.Data = f(x, ...)), attributes(x)))
.lapply <- .preserve.attributes(lapply)
.as.data.frame <- .preserve.attributes(as.data.frame)

# set the row/column values from the file (Origin assumes uniform grid)
.matrix <- function(x) {
	d <- attr(x, 'dimensions')
	dimnames(x) <- list(
		seq(d[4], d[2], length.out = nrow(x)),
		seq(d[3], d[1], length.out = ncol(x))
	)
	t(x)
}

# Long name is easy because it comes first in the \r\n-separated string
# of Long Name, Comment and Units (or whatever).
.get.long.name <- function(lst, name) {
	ret <- sub('\r\n.*$', '', comment(lst)[names(lst) == name])
	if (ret == "") name else ret
}

.expand.tree <- function(tree, lst) lapply(
	setNames( # rename everything to long names, if possible
		tree,
		vapply(
			names(tree),
			# folders (lists here) already have their long names, while
			# other objects must be explicitly renamed if long name exists
			function(n) if (is.list(tree[[n]])) n else .get.long.name(lst, n),
			""
		)
	),
	function(x) switch(typeof(x),
		list = .expand.tree(x, lst),
		character = lst[[x]]
	)
)

# Get the attr-annotated list-of-whatever from the C++ code and transform it
# into usual R data structures.
read.opj <- function(file, encoding = 'latin1', tree = FALSE, ...) {
	ret <- .lapply(
		read_opj(file, encoding, tree), function(x) if (is.list(x)) switch(
			attr(x, 'type'),
			spreadsheet = .as.data.frame(x, ...),
			matrix = .lapply(x, .matrix),
			excel = .lapply(x, .as.data.frame, ...),
		) else x
	)
	if (tree) .expand.tree(attr(ret, 'tree'), ret) else ret
}
