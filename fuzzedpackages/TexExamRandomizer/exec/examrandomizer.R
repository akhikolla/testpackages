#!/usr/local/bin/Rscript
#### TexExamRandomizer SCRIPT WITH JSONPARSER ####
# For use with exams
#### Generating Exams
# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017
#
#
#
# Options:
# Look at ?jsonexamparser for the full list of options.
#
#
library(assertthat)
library(TexExamRandomizer)
library(optparse)

#### RUNNING CODE ####
#

option_list <- list(
    make_option(
        c("--file"),
        action = "store",
        default = NULL,
        type = 'character',
        help = "Filename of the Tex File"
    ),
    make_option(
        c("--table"),
        action = "store",
        default = NULL,
        type = 'character',
        help = "Filename of the table to break down. It overwrites the values written on the file"
    ),
    make_option(
        c("-n", "--noutput"),
        action = "store",
        default = NULL,
        type = "integer",
        help = "Number of output Versions"
    ),
    make_option(
        c("-q", "--nquestions"),
        action = "store",
        default = NULL,
        type = "character",
        help = "Number of output questions"
    ),
    make_option(
        c("-s", "--seed"),
        action = "store",
        default = NULL,
        type = "integer",
        help = "Seed for any randomization done"
    ),
    make_option(
        c("--homework"),
        action = "store_true",
        default = F,
        type = "logical",
        help = "If this option is given. It uses GenerateHomework, which only applies the personalization options, it doesn't try to randomize the order of the document"
    ),
    make_option(
        c("-c", "--compile"),
        action = "store_true",
        default = F,
        type = "logical",
        help = "Should the output folder be compiled or not"
    ),
    make_option(
        c("--xelatex"),
        action = "store_true",
        default = F,
        type = "logical",
        help = "Should we use xelatex to compile or not"
    ),
    make_option(
        c("-d", "--debug"),
        action = "store_true",
        default = F,
        type = "logical",
        help = "If debugging, it doesn't remove auxiliary files"
    )
)


#### PARSING OPTIONS ####
####
opt <-
    parse_args(
        OptionParser(option_list = option_list),
        positional_arguments = T
    )

if (opt$options$homework) {
    jsonhwparser(opt)
} else {
    jsonexamparser(opt)
}






