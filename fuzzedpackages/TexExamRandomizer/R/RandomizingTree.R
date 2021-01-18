# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017
#### Randomizing layers
####
### NOTE: Right now the naming convention is hardcoded into the regex's, Is there any reason to change that?

#' @title GetLayerSampleIndexes
#' @description Function to randomize the names of a list to sample, as provided by \code{\link{StructureDocument}}. It takes the names of the list, not the list itself, and it provides the indexes needed to resample the original list judging by those names. It doesn't output a resampled list.
#'
#'
#' @details Following the prescription from \code{\link{StructureDocument}}, it keeps the "prior_to" and "post_to" parts fixed. And within each section it keeps the "begin_" and "end_" parts fixed. Then, it resamples everything within each section, and afterwards resamples the order of the sections if \code{sampleSectionOrder} is true
#'
#'  Note how the names of the list are expected to represent the structure described in \code{\link{StructureDocument}}
#
#' @param NamesOfListToSample
#'     The names of each element on the list we want to sample, with the original order. The usual use would be \code{NamesOfListToSample = names(list)}.
#'
#'     Being specific, the names that are provided will be matched with the regular expression \code{"^[[:digit:]]+_.+[[:digit:]]+_"}.
#' @param sampleSectionOrder
#'     Should it also move around the sections or not? Look at the details for a more detailed explanation
#' @param isRandomized  If this is set to false, it will not randomize the list, it will just output \code{1:length(OfNamesOfListToSample)}. It handles properly when the list is of length 0, by outputting \code{integer(0)}
#
#
#' @return
#' An integer vector, that would provide a random reordering of the list.
#'
#' @keywords internal
GetLayerSampleIndexes <- function(section_vector, command_vector, sampleSectionOrder = FALSE, randomizeSection = TRUE) {

    # assertthat::assert_that(length(section_vector) == length(command_vector)) This is assumed to be true in here

    IsCorrectlySetUp <- function(toSample) {
        #should be a solid chunk, as it should always leave things on both sides
        return(
            assertthat::are_equal(
                seq.int(from = toSample[1], to = utils::tail(toSample, 1)),
                toSample
            )
        )


    }

    SampleOneSection <- function(cmd_numbers_of_section, line_numbers_of_section, randomizeSection = TRUE) {

        if (!randomizeSection) {
            return(line_numbers_of_section)
        }



        lines_to_sample <- line_numbers_of_section[cmd_numbers_of_section != 0L]

        if (length(lines_to_sample) == 0L) {
            warning("Layer without sublevels as expected on the function call")
            return(line_numbers_of_section)
        } else if (length(lines_to_sample) == 1L) {
            # Only one element, there is nothign to sample
            return(line_numbers_of_section)
        }



        assertthat::assert_that(IsCorrectlySetUp(lines_to_sample))

        start_sample_line <- lines_to_sample[1L]
        end_sample_line <- utils::tail(lines_to_sample, 1L)

        fixed_lines_before <-
            line_numbers_of_section[cmd_numbers_of_section == 0L & line_numbers_of_section < start_sample_line]

        if (end_sample_line == utils::tail(line_numbers_of_section, 1)) {
            fixed_lines_after <- integer(0L)
        } else {
            fixed_lines_after <-
                line_numbers_of_section[cmd_numbers_of_section == 0L & line_numbers_of_section > end_sample_line]
        }



        #CAREFULL WITH SAMPLE. sample(23) is really sample(1:23), while sample (23:24) is simply sample(23:24)... so it will create unintentional behaviour if one doesn't take care (like me)
        #
        return(
            c(
                fixed_lines_before,
                sample(start_sample_line:end_sample_line),
                fixed_lines_after
            )
        )

    }

    # if (!isRandomized) {
    #     return(seq_along(section_vector))
    # }


    nLines <- length(section_vector)

    #   which Lines are part of our structure 1_questions_1_question or so on
    LinesInSections <-
        which(
            section_vector != 0
        )

    if (length(LinesInSections) == 0) {
        # No inside sections
        return(seq_along(section_vector))
    }
    assertthat::assert_that(IsCorrectlySetUp(LinesInSections))
    #There shouldn't be anything inbetween the LinesInSections




    beforeSample <- seq_len(LinesInSections[1] - 1)

    lastLineInSections <- utils::tail(LinesInSections, 1)

    if (lastLineInSections < nLines) {
        # Whatever it is after we start having all the section with our structure
        afterSample <-  seq.int(from = lastLineInSections + 1, to = nLines)
    } else {
        afterSample <- integer()
    }

    uniqueSectionNumber <- unique(section_vector[section_vector != 0])

    #Removes those matches that didn't return a number
    #
    if (sampleSectionOrder && length(uniqueSectionNumber) != 1) {
        uniqueSectionNumber <- sample(uniqueSectionNumber)
    }


    outputSample <- integer()

    for (this_section in uniqueSectionNumber) {
        linesofThisSection <-
            which(
                section_vector == this_section
            )

        outputSample <-
            c(
                outputSample,
                SampleOneSection(
                    command_vector[linesofThisSection],
                    linesofThisSection,
                    randomizeSection
                )
            )

    }

    return(c(beforeSample, outputSample, afterSample))

}



#' @title Randomizing documents.
#' @description Function to randomize a Document, as created by \code{\link{StructureDocument}}.
#'
#' It Randomizes each layer according to the prescriptions involved in the internal function \code{\link{GetLayerSampleIndexes}}. Which, in summary, randomizes each section inside, and then randomizes the orders of the sections.
#'
#' \strong{Important note:} One must provide to this function the \emph{document} part of the structure. Since \code{\link{StructureDocument}} provides as the outer most layer a split between the preamble and the document, one must just supply the document part to this function, (or a subsection of it).
#'
#' @details It keeps randomizing recursively inner layers of the structure until it runs out of elements on the logical vectors \code{isSectionReordered.vector} and \code{isLayerRandomized.vector}.
#'
#' A "section" denotes the content within a begin-end environment in the document.
#' Each section is then assumed to be divided in a beginning and end parts, that should be fixed in place, and the parts denoted by the command \code{\\cmdName} as explained on \code{\link{StructureDocument}}.
#'
#' We will denote those parts as "items."
#' Analogously to itemize environments in 'LaTeX'.
#'
#' The purpose of this function is therefore to randomize the items from the structure, fixing the begin and end parts within a section. And then to reorder each section while keeping the pre- and post- parts fixed, and to do so recursively until we exhaust the \code{isLayerRandomized.vector}
#'
#' \code{isSectionReordered.vector} specifies whether to order sections for a certain depth, while \code{isLayerRandomized.vector} specifies whether to order the items within a section of that same depth.
#'
#' In some cases you may want to reorder the sections, for example, using the examdesign class. Over there, questions use the begin-end question format.
#'
#' In others cases you may want to preserve the order of sections while still modifying the order of the items, like when you are using the exam class, or when creating your own list of questions with an \code{\\itemize} environment.
#'
#' For efficiency, if you don't want to randomize to the full depth of your tree, just make those logical vectors of your desired length, rather than making them of length \eqn{n} and then setting every layer after the last one you want to randomize to false. That will prevent the program from walking down the whole tree checking everything.
#'
#' @param Document Document to randomize, as generated by \code{\link{StructureDocument}}. The names of the structure are used to determine how to randomize the document.
#' @param isSectionReordered.vector Logical vector, specifying if the order of sections should be also randomized at a certain depth level.
#'
#' \strong{Note} that if \code{isLayerRandomized} is set to false for a certain layer, \code{isSectionReorder} will have no effect. (Probably this isn't the best behaviour)
#' @param isLayerRandomized.vector Logical vector, specifying if you should randomize the order of the items, (denoted by \code{\\cmdName}) or not at a certain depth level.
#'
#' This vector should have the same length as the depth at most, otherwise it will raise an error if you try to "dig deeper than it can". And \emph{isSectionReordered.vector} should have matching elements for each element of \emph{isLayerRandomized.vector}
#'
#' (Maybe we could change this to a warning instead? To allow for structures with different depths within different branches of the tree)
#'
#'
#' @return A document structure, as provided by \code{\link{StructureDocument}}.
#'
#' However, the names of the structure will no longer be sequential,
#' the naming convention in the new structure will refer to the original structure that was inputted into this function. Which is very useful when you want to keep track of where things have moved.
#' @seealso \code{\link{StructureDocument}}, TODO: Add reference to extracting info functions
#' @export
#' @examples
#'
#' rndDoc <- RandomizeDocument(
#'     TexExamRandomizer::testdoc$document,
#'     c(FALSE,TRUE),
#'     c(TRUE, TRUE)
#'     )

RandomizeDocument <- function(Document, isSectionReordered.vector, isLayerRandomized.vector) {


    #It should mantain 0s where it is, and then it changes consequtive numbers into further and further moves

    reformat <- function(vector) {

        if (length(vector) == 0) {
            return(NULL)
        }
        reformatted_vector <- integer(length(vector))
        value_index <- 1L
        last_value <- 0L

        for (i in seq_along(vector)) {
            if (vector[i] == 0L) {
                reformatted_vector[i] <- 0L
            } else if (vector[i] != last_value) {
                last_value <- last_value + 1
                reformatted_vector[i] <- last_value
            } else {
                reformatted_vector[i] <- last_value
            }

        }

        return(reformatted_vector)

    }


    #This function does the recursive part of the algorithm.
    # Relies on     eval.parent(substitute( <- )) to imitate a "passing by reference" behaviour,
    # since this function updates ListWithLayers to the new desired one.
    # Using a proper pass by reference could probably optimize this function... if there is any.
    RecurviseRandomization <- function(ListWithLayers, left.isSectionReordered, left.layers) {


        #If this assertion isn't true, that means we are trying to dig deeper into the structure than we have actually structured it to.
        assertthat::assert_that(
            is.list(ListWithLayers),
            msg = "ListWithLayers wasn't a list.\nAre you sure you are not trying to traverse this document further down than it has been structured to?"
        )


        temp <- ListWithLayers

        section_before <- attr(ListWithLayers, "section")
        command_before <- attr(ListWithLayers, "command")
        assertthat::assert_that(
            !is.null(section_before),
            !is.null(command_before),
            msg = "The list requires a section and command attribute"
        )
        assertthat::assert_that(
            length(ListWithLayers) == length(section_before),
            length(section_before) == length(command_before)
            )
        newIndexes <- GetLayerSampleIndexes(section_before, command_before, left.isSectionReordered[1], left.layers[1])

        if (length(left.layers) > 1) {
            for (i in seq_along(ListWithLayers)) {
                if (command_before[i] != 0) {
                    RecurviseRandomization(temp[[i]], left.isSectionReordered[-1],  left.layers[-1])
                }
            }
        }
        section_original <- attr(ListWithLayers, "section")[newIndexes]
        command_original <- attr(ListWithLayers, "command")[newIndexes]
        section <- reformat(section_original)
        command <- reformat(command_original)

        eval.parent(
            substitute(
                ListWithLayers <-
                    structure(
                        temp[newIndexes],
                        section = section,
                        command = command,
                        section_original = section_original,
                        command_original = command_original
                        )
                )
            )
    }

    # Actual code of the function
    #
    assertthat::assert_that(length(isSectionReordered.vector) >= length(isLayerRandomized.vector))

    RecurviseRandomization(Document, isSectionReordered.vector, isLayerRandomized.vector)

    return(Document)


}
# If you want to remove questions, do it after you have randomized it, just pick the last 5 and get rid of them or something like that, that will make everything a lot simpler.
