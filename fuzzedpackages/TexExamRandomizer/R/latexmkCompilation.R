# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017
#### LateX Folder Compilation Codes ####
####
####

#' @title Compiling function
#'
#' @description This function calls latexmk, which must be part of the system commands, a directory where tex files are found and outputs their pdf and other things in the pdf.dir.out
#' The functions \code{\link{CompileLatexDirEXAM}} and \code{\link{CompileLatexDirHW}} are identical wrappers of the same function, \code{\link{CompileLatexDir}}. Do not use them, they are just kept for "backwards" compatibility
#'
#' @details Write the tex files relative paths to other files as to be read from the directory in which latex.dir.in is found
#' This function is intended to be use to compile a bunch of files which are stemmed from an original one. That is why the directory
#'
#' @param pdf.dir.out Directory where the pdf output will be sent to
#' @param latex.dir.in Directory where all the tex files are found.
#' @param engine: Engine to use when compiling. Currently the options are \code{xelatex}, \code{lualatex}, \code{latex} and \code{pdflatex}
#'
#'        \code{xelatex} is the default value. However, if the value is not recognized, \code{pdflatex} is used instead.
#
#' @param compile.dir: Directory from which compilation is invoked, if not specified, the directory we are compiling will be from where we do it. (This is specially usefull since we want to mantain the same relative paths from the main file).
#' @return None
#' @author Alejandro Recuenco \email{alejandrogonzalezrecuenco@@gmail.com}
#' @export
#'
#' @family Compilation functions
#' @examples
#'
#' input_folder <- system.file(
#'     "extdata",
#'     "ExampleTexDocuments",
#'     package = "TexExamRandomizer")
#'
#'
#' TexExamRandomizer::CompileLatexDir(
#'     pdf.dir.out = tempdir(),
#'     engine= "pdf",
#'     latex.dir.in = input_folder,
#'     extracmdoptions = "-time")
#'

CompileLatexDir <- function(pdf.dir.out, latex.dir.in, engine = "xelatex", compile.dir = NULL, extracmdoptions = NULL){

  full_pdf_output_dir <- normalizePath(pdf.dir.out)
  full_tex_input_dir <- normalizePath(latex.dir.in)


  if        (engine == "xelatex")  {
      cmd_option   <- "-xelatex"
  } else if (engine == "lualatex") {
      cmd_option   <- "-lualatex"
  } else if (engine == "latex")    {
      cmd_option   <- "-latex"
  } else                           {
      cmd_option   <- "-pdf"
  }

  # TODO: try to incorporate match.arg to simplify the code here

  cmdarg_outputdir <- sprintf("-outdir=%s", shQuote(full_pdf_output_dir))

  cmdarg_input_tex_files <- shQuote(
      list.files(path = full_tex_input_dir, pattern = "\\.tex$", full.names = T)
  )

  extracmdoptions <- c(cmd_option, extracmdoptions)
  othercmdoption <- "-quiet"


  cmdArgs <- c(cmdarg_outputdir,
               extracmdoptions,
               othercmdoption,
               cmdarg_input_tex_files
               )
  if (is.null(compile.dir)) {
    compile.dir <- latex.dir.in
  }

  owd <- getwd() #Original Working Directory
  if (!is.null(owd)) {
    on.exit(setwd(owd), add = TRUE)
    setwd(compile.dir)
    #After the on.exit, if it causes an error setwd(owd) is still called, do not place the setwd before the on.exit.
  } else {
    warning("Current working directory is unknown, can't change the compilation directory and return to the same directory afterwards...")
  }




  system2(
    command = "latexmk",
    args = cmdArgs,
    wait = TRUE
  );
}



#' @rdname  CompileLatexDir
#' @keywords internal
#' @family Compilation functions
CompileLatexDirEXAM <- function(pdf.dir.out, latex.dir.in, engine = "xelatex", compile.dir = NULL, extracmdoptions = NULL){
  CompileLatexDir(pdf.dir.out, latex.dir.in, engine = engine, compile.dir = compile.dir, extracmdoptions = extracmdoptions)
}

#' @rdname CompileLatexDir
#' @keywords internal
#' @family Compilation functions
CompileLatexDirHW <- function(pdf.dir.out, latex.dir.in, engine = "xelatex", compile.dir = NULL, extracmdoptions = NULL){
    CompileLatexDir(pdf.dir.out, latex.dir.in, engine = engine, compile.dir = compile.dir, extracmdoptions = extracmdoptions)
}
