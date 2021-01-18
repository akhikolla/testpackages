
#' chunker class
#' @name chunker-class
#' @slot pointer externalptr object
#' @keywords internal

setClass( "chunker", representation( pointer = "externalptr", attr = "list"))