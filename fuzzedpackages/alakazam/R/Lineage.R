# Ig lineage reconstruction via maximum parsimony

#' @include Classes.R
NULL


#### Preprocessing functions ####

#' Generate a ChangeoClone object for lineage construction
#' 
#' \code{makeChangeoClone} takes a data.frame with AIRR or Change-O style columns as input and 
#' masks gap positions, masks ragged ends, removes duplicates sequences, and merges 
#' annotations associated with duplicate sequences. It returns a \code{ChangeoClone} 
#' object which serves as input for lineage reconstruction.
#' 
#' @param    data         data.frame containing the AIRR or Change-O data for a clone. See Details
#'                        for the list of required columns and their default values.
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing observed DNA sequences. All 
#'                        sequences in this column must be multiple aligned.
#' @param    germ         name of the column containing germline DNA sequences. All entries 
#'                        in this column should be identical for any given clone, and they
#'                        must be multiple aligned with the data in the \code{seq} column.
#' @param    v_call       name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    j_call        name of the column containing J-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a 
#'                        numeric value. All entries in this column should be identical 
#'                        for any given clone.
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    mask_char    character to use for masking and padding.
#' @param    max_mask     maximum number of characters to mask at the leading and trailing
#'                        sequence ends. If \code{NULL} then the upper masking bound will 
#'                        be automatically determined from the maximum number of observed 
#'                        leading or trailing Ns amongst all sequences. If set to \code{0} 
#'                        (default) then masking will not be performed.
#' @param    pad_end      if \code{TRUE} pad the end of each sequence with \code{mask_char}
#'                        to make every sequence the same length.
#' @param    text_fields  text annotation columns to retain and merge during duplicate removal.
#' @param    num_fields   numeric annotation columns to retain and sum during duplicate removal.
#' @param    seq_fields   sequence annotation columns to retain and collapse during duplicate 
#'                        removal. Note, this is distinct from the \code{seq} and \code{germ} 
#'                        arguments, which contain the primary sequence data for the clone
#'                        and should not be repeated in this argument.
#' @param    add_count    if \code{TRUE} add an additional annotation column called 
#'                        \code{collapse_count} during duplicate removal that indicates the 
#'                        number of sequences that were collapsed.
#' @param    verbose      passed on to \code{collapseDuplicates}. If \code{TRUE}, report the 
#'                        numbers of input, discarded and output sequences; otherwise, process
#'                        sequences silently.                        
#'
#' @return   A \link{ChangeoClone} object containing the modified clone.
#'
#' @details
#' The input data.frame (\code{data}) must columns for each of the required column name 
#' arguments: \code{id}, \code{seq}, \code{germ}, \code{v_call}, \code{j_call}, 
#' \code{junc_len}, and \code{clone}.  The default values are as follows:
#' \itemize{
#'   \item  \code{id       = "sequence_id"}:         unique sequence identifier.
#'   \item  \code{seq      = "sequence_alignment"}:  IMGT-gapped sample sequence.
#'   \item  \code{germ     = "germline_alignment"}:  IMGT-gapped germline sequence.
#'   \item  \code{v_call   = "v_call"}:              V segment allele call.
#'   \item  \code{j_call   = "j_call"}:              J segment allele call.
#'   \item  \code{junc_len = "junction_length"}:     junction sequence length.
#'   \item  \code{clone    = "clone_id"}:            clone identifier.
#' }
#' Additional annotation columns specified in the \code{text_fields}, \code{num_fields} 
#' or \code{seq_fields} arguments will be retained in the \code{data} slot of the return 
#' object, but are not required. If the input data.frame \code{data} already contains a 
#' column named \code{sequence}, which is not used as the \code{seq} argument, then that 
#' column will not be retained.
#' 
#' The default columns are IMGT-gapped sequence columns, but this is not a requirement. 
#' However, all sequences (both observed and germline) must be multiple aligned using
#' some scheme for both proper duplicate removal and lineage reconstruction. 
#'
#' The value for the germline sequence, V-segment gene call, J-segment gene call, 
#' junction length, and clone identifier are determined from the first entry in the 
#' \code{germ}, \code{v_call}, \code{j_call}, \code{junc_len} and \code{clone} columns, 
#' respectively. For any given clone, each value in these columns should be identical.
#'  
#' @seealso  Executes in order \link{maskSeqGaps}, \link{maskSeqEnds}, 
#'           \link{padSeqEnds}, and \link{collapseDuplicates}. 
#'           Returns a \link{ChangeoClone} object which serves as input to
#'           \link{buildPhylipLineage}.
#' 
#' @examples
#' # Example data
#' db <- data.frame(sequence_id=LETTERS[1:4],
#'                  sequence_alignment=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
#'                  germline_alignment="CCCCAGGG",
#'                  v_call="Homsap IGKV1-39*01 F",
#'                  j_call="Homsap IGKJ5*01 F",
#'                  junction_length=2,
#'                  clone_id=1,
#'                  c_call=c("IGHM", "IGHG", "IGHG", "IGHA"),
#'                  duplicate_count=1:4,
#'                  stringsAsFactors=FALSE)
#' 
#' 
#'  # Without end masking
#'  makeChangeoClone(db, text_fields="c_call", num_fields="duplicate_count")
#' 
#'  # With end masking
#'  makeChangeoClone(db, max_mask=3, text_fields="c_call", num_fields="duplicate_count")
#'
#' @export
makeChangeoClone <- function(data, id="sequence_id", seq="sequence_alignment", 
                             germ="germline_alignment", v_call="v_call", j_call="j_call",
                             junc_len="junction_length", clone="clone_id", mask_char="N",
                             max_mask=0, pad_end=FALSE, text_fields=NULL, num_fields=NULL, seq_fields=NULL,
                             add_count=TRUE, verbose=FALSE) {
    # Check for valid fields
    check <- checkColumns(data, c(id, seq, germ, v_call, j_call, junc_len, clone, 
                                  text_fields, num_fields, seq_fields))
    if (check != TRUE) { stop(check) }
    
    # Replace gaps with Ns and masked ragged ends
    tmp_df <- data[, c(id, seq, text_fields, num_fields, seq_fields)]
    tmp_df[[seq]] <- maskSeqGaps(tmp_df[[seq]], mask_char=mask_char, outer_only=FALSE)
    tmp_df[[seq]] <- maskSeqEnds(tmp_df[[seq]], mask_char=mask_char, max_mask=max_mask, trim=FALSE)
    
    # Pad ends
    if (pad_end) {
        tmp_df[[seq]] <- padSeqEnds(tmp_df[[seq]], pad_char=mask_char)
    }
    
    seq_len <- stringi::stri_length(tmp_df[[seq]])
    if (any(seq_len != seq_len[1])) {
        len_message <- paste0("All sequences are not the same length for data with first ", 
                              id, " = ", tmp_df[[id]][1], ".")
        if (!pad_end) {
            len_message <- paste(len_message, 
                                 "Consider specifying pad_end=TRUE and verify the multiple alignment.")
        } else {
            len_message <- paste(len_message,
                                 "Verify that all sequences are properly multiple-aligned.")
        }
        stop(len_message)
    }
    
    # Remove duplicates
    tmp_df <- collapseDuplicates(tmp_df, id=id, seq=seq, text_fields=text_fields, 
                                 num_fields=num_fields, seq_fields=seq_fields,
                                 add_count=add_count, verbose=verbose)
    
    # Define return object
    tmp_names <- names(tmp_df)
    if ("sequence" %in% tmp_names & seq != "sequence") {
        tmp_df <- tmp_df[, tmp_names != "sequence"]
        tmp_names <- names(tmp_df)
    }
    names(tmp_df)[tmp_names == seq] <- "sequence"
    names(tmp_df)[tmp_names == id] <- "sequence_id"
    clone <- new("ChangeoClone", 
                 data=as.data.frame(tmp_df),
                 clone=as.character(data[[clone]][1]),
                 germline=maskSeqGaps(data[[germ]][1], mask_char=mask_char, outer_only=FALSE), 
                 v_gene=getGene(data[[v_call]][1]), 
                 j_gene=getGene(data[[j_call]][1]), 
                 junc_len=data[[junc_len]][1])
    
    return(clone)
}


#### PHYLIP functions ####

# Create PHYLIP input files in a temporary folder
#
# @param   clone  a ChangeoClone object
# @param   path   a directory to store the write the output files to
# @return  a named vector translating sequence_id (names) to PHYLIP taxa (values)
writePhylipInput <- function(clone, path) {
    # Define PHYLIP columns
    nseq <- nrow(clone@data)
    v1 <- c(sprintf('%-9s', nseq + 1),
            sprintf("%-9s", "Germline"), 
            sprintf("SAM%-6s", 1:nseq))
    v2 <- c(stringi::stri_length(clone@germline),
            clone@germline, 
            clone@data[["sequence"]])
    phy_df <- data.frame(v1, v2, stringsAsFactors=F)
    
    # Define names vector mapping taxa names to original sequence identifiers
    id_map <- setNames(gsub("^\\s+|\\s+$", "", v1[-(1:2)]), clone@data[["sequence_id"]])
    
    # Create PHYLIP input file
    write.table(phy_df, file=file.path(path, "infile"), 
                quote=F, sep=" ", col.names=F, row.names=F)    
    
    return(id_map)
}


# Run PHYLIP dnapars application
#
# @param   path          temporary directory containing infile.
# @param   phylip_exec  path to the dnapars executable.
# @param   verbose       if TRUE suppress phylip console output.
# @return  TRUE if phylip ran successfully and FALSE otherwise
runPhylip <- function(path, phylip_exec, verbose=FALSE) {
    # Expand shell variables
    phylip_exec <- path.expand(phylip_exec)
    
    # Remove old files
    if (file.exists(file.path(path, "outfile"))) { file.remove(file.path(path, "outfile")) }
    if (file.exists(file.path(path, "outtree"))) { file.remove(file.path(path, "outtree")) }    
    
    # Set platform specific options
    if (.Platform$OS.type == "windows") { 
        quiet_params <- list(ignore.stdout=TRUE, ignore.stderr=TRUE)
        invoke <- shell
    } else { 
        quiet_params <- list(stdout=FALSE, stderr=FALSE)
        invoke <- system2
    } 
    
    # Set dnapars options
    phy_options <- c("S", "Y", "I", "4", "5", ".")
    params <- list(phylip_exec, input=c(phy_options, "Y"), wait=TRUE)
    if (!verbose) {
        params <- append(params, quiet_params)
    }
    
    # Call phylip
    wd <- getwd()
    setwd(path)
    status <- tryCatch(do.call(invoke, params), error=function(e) e)
    setwd(wd)
    
    # Return TRUE if phylip ran successfully
    invisible(status == 0)
}


# Reads in the PHYLIP outfile
#
# @param   path  the temporary folder containing the dnapars outfile
# @return  a character vector with each item as a line in the outfile
readPhylipOutput <- function(path) {
    phylip_out <- scan(file.path(path, "outfile"), what="character", sep="\n", 
                       blank.lines.skip=FALSE, strip.white=FALSE, quiet=TRUE)
    return(phylip_out)
}


# Test for successful PHYLIP dnapars run by checking the outfile
#
# @param   phylip_out  a character vector returned by readPhylipOut
# @return  TRUE if trees built 
#          FALSE if no trees built
checkPhylipOutput <- function(phylip_out) {
    # Check for failed tree build
    result <- !(any(grepl('-1 trees in all found', phylip_out)))
    
    return(result)
}


# Extracts inferred sequences from PHYLIP dnapars outfile
#
# @param   phylip_out   a character vector returned by readPhylipOutput
# @return  a list containing an id vector, a sequence vector and an annotation data.frame
getPhylipInferred <- function(phylip_out) {
    # Process dnapars output
    seq_start <- min(grep("From\\s+To\\s+Any Steps\\?\\s+State at upper node", 
                          phylip_out, perl=T, fixed=F))
    seq_empty <- grep("^\\s*$", phylip_out[seq_start:length(phylip_out)], perl=T, fixed=F)
    seq_len <- seq_empty[min(which(seq_empty[-1] == (seq_empty[-length(seq_empty)] + 1)))]
    seq_block <- paste(phylip_out[(seq_start + 2):(seq_start + seq_len - 2)], collapse="\n")
    seq_df <- read.table(textConnection(seq_block), as.is=T, fill=T, blank.lines.skip=F)
    
    # Correct first line of block and remove blank rows
    fix.row <- c(1, which(is.na(seq_df[,1])) + 1)
    end_col <-  ncol(seq_df) - 2
    #seq_df[fix.row, ] <- cbind(0, seq_df[fix.row, 1], "no", seq_df[fix.row, 2:5], stringsAsFactors=F)
    seq_df[fix.row, ] <- data.frame(cbind(0, seq_df[fix.row, 1], "no", seq_df[fix.row, 2:end_col]), stringsAsFactors=F)
    if (length(fix.row)>1) {
        seq_df <- seq_df[-(fix.row[-1] - 1), ]
    }
    
    # Create data.frame of inferred sequences
    inferred_num <- unique(grep("^[0-9]+$", seq_df[, 2], value=T))
    inferred_seq <- sapply(inferred_num, function(n) { paste(t(as.matrix(seq_df[seq_df[, 2] == n, -c(1:3)])), collapse="") })
    
    if (length(inferred_num)>0) {
        return(data.frame(sequence_id=paste0("Inferred", inferred_num), 
                          sequence=inferred_seq, stringsAsFactors = FALSE))
    }
    data.frame(sequence_id=c(), sequence=c(), stringsAsFactors = FALSE)
}


# Extracts graph edge list from a PHYLIP dnapars outfile
#
# @param   phylip_out  character vector returned by readPhylipOutput
# @param   id_map      named vector of PHYLIP taxa names (values) to sequence 
#                      identifiers (names) that will be translated. If NULL
#                      no taxa name translation is performed
# @return  a data.frame of edges with columns (from, to, weight)
getPhylipEdges <- function(phylip_out, id_map=NULL) {
    # Process dnapars output
    edge_start <- min(grep('between\\s+and\\s+length', phylip_out, 
                           perl=TRUE, fixed=FALSE))
    edge_len <- min(grep('^\\s*$', phylip_out[edge_start:length(phylip_out)], 
                         perl=TRUE, fixed=FALSE))
    edge_block <- paste(phylip_out[(edge_start + 2):(edge_start + edge_len - 2)], collapse='\n')
    edge_df <- read.table(textConnection(edge_block), col.names=c('from', 'to', 'weight'), 
                          as.is=TRUE)

    # Modify inferred taxa names to include "Inferred"
    inf_map <- unique(grep("^[0-9]+$", c(edge_df$from, edge_df$to), value=T))
    names(inf_map) <- paste0("Inferred", inf_map)
    edge_df$from <- translateStrings(edge_df$from, inf_map)
    edge_df$to <- translateStrings(edge_df$to, inf_map)
    
    if (!is.null(id_map)) {
        # Reassign PHYLIP taxa names to sequence IDs
        edge_df$from <- translateStrings(edge_df$from, id_map)
        edge_df$to <- translateStrings(edge_df$to, id_map)
    }
    
    return(edge_df)
}


# Modify edges of phylip output
#
# @param   edges     data.frame of edges returned by getPhylipEdges
# @param   clone     a ChangeoClone object containg sequence data
# @param   dist_mat  DNA character distance matrix
# @return  a list of modified edges data.frame and clone object
modifyPhylipEdges <- function(edges, clone, dist_mat=getDNAMatrix(gap=0)) {
    # Move germline to root position
    germ_idx <- which(edges$to == "Germline")
    edges[germ_idx, c('from', 'to')] <- edges[germ_idx, c('to', 'from')]

    # Calculate edge mutations
    for (i in 1:nrow(edges)) {
        if (edges$from[i] == "Germline") {
            seq1 <- clone@germline
        } else {
            seq1 <- clone@data[["sequence"]][clone@data[["sequence_id"]] == edges$from[i]]
        }
        seq2 <- clone@data[["sequence"]][clone@data[["sequence_id"]] == edges$to[i]]
        edges$weight[i] <- seqDist(seq1, seq2, dist_mat)        
    }
    
    # Find rows zero weight edges with inferred parent nodes
    remove_row <- which(edges$weight == 0 & 
                        edges$from != "Germline" & 
                        grepl('^Inferred\\d+$', edges$from))
    
    # Replace inferred parent nodes with child nodes when edge weight is zero
    while (length(remove_row) > 0) {
        # Remove first node with zero distance to parent
        r <- remove_row[1]
        r_idx <- which(edges[c('from', 'to')] == edges$from[r], arr.ind=T)
        edges[r_idx] <- edges$to[r]
        
        # Recalculate edge weights for modified rows
        r_mod <- r_idx[, 1][r_idx[, 1] != r]
        for (i in r_mod) {
            if (edges$from[i] == "Germline") {
                seq1 <- clone@germline
            } else {
                seq1 <- clone@data[["sequence"]][clone@data[["sequence_id"]] == edges$from[i]]
            }
            seq2 <- clone@data[["sequence"]][clone@data[["sequence_id"]] == edges$to[i]]
            edges$weight[i] <- seqDist(seq1, seq2, dist_mat)      
        }
        
        # Remove row
        edges <- edges[-r, ]
        
        # Re-determine rows to remove
        remove_row <- which(edges$weight == 0 & 
                            edges$from != "Germline" & 
                            grepl('^Inferred\\d+$', edges$from))      
    }
    
    # Remove rows from clone
    keep_clone <- clone@data[["sequence_id"]] %in% unique(c(edges$from, edges$to))
    clone@data <- as.data.frame(clone@data[keep_clone, ])
    
    return(list(edges=edges, clone=clone))
}

# Convert edge data.frame and clone object to igraph graph object
#
# @param   edges  data.frame of edges returned by getPhylipEdges
# @param   clone  a ChangeoClone object containg sequence data
# @return  an igraph graph object
phylipToGraph <- function(edges, clone) {
    # Create igraph object
    g <- igraph::graph_from_data_frame(edges, directed=T)
    
    # Add germline sequence
    germ_idx <- which(igraph::V(g)$name == "Germline")
    g <- igraph::set_vertex_attr(g, "sequence", index=germ_idx, clone@germline)
    
    # Add sample sequences and names
    clone_idx <- match(clone@data[["sequence_id"]], igraph::V(g)$name) 
    g <- igraph::set_vertex_attr(g, "sequence", index=clone_idx, clone@data[["sequence"]])
    
    # Add annotations
    ann_fields <- names(clone@data)[!(names(clone@data) %in% c("sequence_id", "sequence"))]
    for (n in ann_fields) {
        g <- igraph::set_vertex_attr(g, n, index=germ_idx, NA)
        g <- igraph::set_vertex_attr(g, n, index=clone_idx, clone@data[[n]])
    }
    
    # Add edge and vertex labels
    igraph::V(g)$label <- igraph::V(g)$name
    igraph::E(g)$label <- igraph::E(g)$weight
    
    # Add graph attributes
    g$clone <- clone@clone
    g$v_gene <- clone@v_gene
    g$j_gene <- clone@j_gene
    g$junc_len <- clone@junc_len
    
    return(g)
}


#' Infer an Ig lineage using PHYLIP
#' 
#' \code{buildPhylipLineage} reconstructs an Ig lineage via maximum parsimony using the 
#' dnapars application of the PHYLIP package.
#' 
#' @param    clone         \link{ChangeoClone} object containing clone data.
#' @param    phylip_exec   absolute path to the PHYLIP dnapars executable.
#' @param    dist_mat      Character distance matrix to use for reassigning edge weights. 
#'                         Defaults to a Hamming distance matrix returned by \link{getDNAMatrix} 
#'                         with \code{gap=0}. If gap characters, \code{c("-", ".")}, are assigned 
#'                         a value of -1 in \code{dist_mat} then contiguous gaps of any run length,
#'                         which are not present in both sequences, will be counted as a 
#'                         distance of 1. Meaning, indels of any length will increase
#'                         the sequence distance by 1. Gap values other than -1 will 
#'                         return a distance that does not consider indels as a special case.
#' @param    rm_temp       if \code{TRUE} delete the temporary directory after running dnapars;
#'                         if \code{FALSE} keep the temporary directory.
#' @param    verbose       if \code{FALSE} suppress the output of dnapars; 
#'                         if \code{TRUE} STDOUT and STDERR of dnapars will be passed to 
#'                         the console.
#'                                                
#' @return   An igraph \code{graph} object defining the Ig lineage tree. Each unique input 
#'           sequence in \code{clone} is a vertex of the tree, with additional vertices being
#'           either the germline (root) sequences or inferred intermediates. The \code{graph} 
#'           object has the following attributes.
#'           
#'           Vertex attributes:
#'           \itemize{
#'             \item  \code{name}:      value in the \code{sequence_id} column of the \code{data} 
#'                                      slot of the input \code{clone} for observed sequences. 
#'                                      The germline (root) vertex is assigned the name 
#'                                      "Germline" and inferred intermediates are assigned
#'                                      names with the format {"Inferred1", "Inferred2", ...}.
#'             \item  \code{sequence}:  value in the \code{sequence} column of the \code{data} 
#'                                      slot of the input \code{clone} for observed sequences.
#'                                      The germline (root) vertex is assigned the sequence
#'                                      in the \code{germline} slot of the input \code{clone}.
#'                                      The sequence of inferred intermediates are extracted
#'                                      from the dnapars output.
#'             \item  \code{label}:     same as the \code{name} attribute.
#'           }
#'           Additionally, each other column in the \code{data} slot of the input 
#'           \code{clone} is added as a vertex attribute with the attribute name set to 
#'           the source column name. For the germline and inferred intermediate vertices,
#'           these additional vertex attributes are all assigned a value of \code{NA}.
#'           
#'           Edge attributes:
#'           \itemize{
#'             \item  \code{weight}:    Hamming distance between the \code{sequence} attributes
#'                                      of the two vertices.
#'             \item  \code{label}:     same as the \code{weight} attribute.
#'           }
#'           Graph attributes:
#'           \itemize{
#'             \item  \code{clone}:     clone identifier from the \code{clone} slot of the
#'                                      input \code{ChangeoClone}.
#'             \item  \code{v_gene}:    V-segment gene call from the \code{v_gene} slot of 
#'                                      the input \code{ChangeoClone}.
#'             \item  \code{j_gene}:    J-segment gene call from the \code{j_gene} slot of 
#'                                      the input \code{ChangeoClone}.
#'             \item  \code{junc_len}:  junction length (nucleotide count) from the 
#'                                      \code{junc_len} slot of the input \code{ChangeoClone}.
#'           }
#'           
#' @details
#' \code{buildPhylipLineage} builds the lineage tree of a set of unique Ig sequences via
#' maximum parsimony through an external call to the dnapars application of the PHYLIP
#' package. dnapars is called with default algorithm options, except for the search option, 
#' which is set to "Rearrange on one best tree". The germline sequence of the clone is used 
#' for the outgroup. 
#' 
#' Following tree construction using dnapars, the dnapars output is modified to allow
#' input sequences to appear as internal nodes of the tree. Intermediate sequences 
#' inferred by dnapars are replaced by children within the tree having a Hamming distance 
#' of zero from their parent node. With the default \code{dist_mat}, the distance calculation 
#' allows IUPAC ambiguous character matches, where an ambiguous character has distance zero 
#' to any character in the set of characters it represents. Distance calculation and movement of 
#' child nodes up the tree is repeated until all parent-child pairs have a distance greater than zero 
#' between them. The germline sequence (outgroup) is moved to the root of the tree and
#' excluded from the node replacement processes, which permits the trunk of the tree to be
#' the only edge with a distance of zero. Edge weights of the resultant tree are assigned 
#' as the distance between each sequence.
#' 
#' @references
#' \enumerate{
#'   \item  Felsenstein J. PHYLIP - Phylogeny Inference Package (Version 3.2). 
#'            Cladistics. 1989 5:164-166.
#'   \item  Stern JNH, Yaari G, Vander Heiden JA, et al. B cells populating the multiple 
#'            sclerosis brain mature in the draining cervical lymph nodes. 
#'            Sci Transl Med. 2014 6(248):248ra107.
#' }
#'   
#' @seealso  Takes as input a \link{ChangeoClone}. 
#'           Temporary directories are created with \link{makeTempDir}.
#'           Distance is calculated using \link{seqDist}. 
#'           See \link{igraph} and \link{igraph.plotting} for working 
#'           with igraph \code{graph} objects. 
#'
#' @examples
#' \dontrun{
#' # Preprocess clone
#' db <- subset(ExampleDb, clone_id == 3138)
#' clone <- makeChangeoClone(db, text_fields=c("sample_id", "c_call"), 
#'                           num_fields="duplicate_count")
#' 
#' # Run PHYLIP and process output
#' phylip_exec <- "~/apps/phylip-3.695/bin/dnapars"
#' graph <- buildPhylipLineage(clone, phylip_exec, rm_temp=TRUE)
#' 
#' # Plot graph with a tree layout
#' library(igraph)
#' plot(graph, layout=layout_as_tree, vertex.label=V(graph)$c_call, 
#'      vertex.size=50, edge.arrow.mode=0, vertex.color="grey80")
#' 
#' # To consider each indel event as a mutation, change the masking character 
#' # and distance matrix
#' clone <- makeChangeoClone(db, text_fields=c("sample_id", "c_call"), 
#'                           num_fields="duplicate_count", mask_char="-")
#' graph <- buildPhylipLineage(clone, phylip_exec, dist_mat=getDNAMatrix(gap=-1), 
#'                             rm_temp=TRUE)
#' }
#' 
#' @export
buildPhylipLineage <- function(clone, phylip_exec, dist_mat=getDNAMatrix(gap=0), 
                               rm_temp=FALSE, verbose=FALSE) {
    # Check clone size
    if (nrow(clone@data) < 2) {
        warning("Clone ", clone@clone, " was skipped as it does not contain at least 
                2 unique sequences")
        return(NULL)
    }
    
    # Check fields
    seq_len = unique(stringi::stri_length(clone@data[["sequence"]]))
    germ_len = ifelse(length(clone@germline) == 0, 0, stringi::stri_length(clone@germline))
    if(germ_len == 0) {
        stop("Clone ", clone@clone, "does not contain a germline sequence.")
    }
    if(length(seq_len) != 1) {
        stop("Clone ", clone@clone, "does not contain sequences of equal length.")
    }
    if(seq_len != germ_len) {
        stop("The germline and input sequences are not the same length for clone ", clone@clone)
    }
    
    # Check dnapars access
    if (file.access(phylip_exec, mode=1) == -1) {
        stop("The file ", phylip_exec, " cannot be executed.")
    }
    
    # Create temporary directory
    temp_path <- makeTempDir(paste0(clone@clone, "-phylip"))
    if (verbose) {
        cat("TEMP_DIR> ", temp_path, "\n", sep="")
    }
    
    # Run PHYLIP
    id_map <- writePhylipInput(clone, temp_path)
    runPhylip(temp_path, phylip_exec, verbose=verbose)
    phylip_out <- readPhylipOutput(temp_path)
    
    # Remove temporary directory
    if (rm_temp) {
        unlink(temp_path, recursive=TRUE)
    }
    
    # Check output for trees
    if (!checkPhylipOutput(phylip_out)) {
        warning('PHYLIP failed to generate trees for clone ', clone)
        return(NULL)
    }
    
    # Extract inferred sequences from PHYLIP output
    inf_df <- getPhylipInferred(phylip_out)
    clone@data <- as.data.frame(dplyr::bind_rows(clone@data, inf_df))

    # Extract edge table from PHYLIP output 
    edges <- getPhylipEdges(phylip_out, id_map=id_map)
    
    # Modify PHYLIP tree to remove 0 distance edges
    mod_list <- modifyPhylipEdges(edges, clone, dist_mat=dist_mat)

    # Convert edges and clone data to igraph graph object
    graph <- phylipToGraph(mod_list$edges, mod_list$clone)
    
    return(graph)
}


#' Convert a tree in ape \code{phylo} format to igraph \code{graph} format.
#' 
#' \code{phyloToGraph} converts a tree in \code{phylo} format to and 
#' \code{graph} format.
#' 
#' @param   phylo  An ape \code{phylo} object.
#' @param   germline  If specified, places specified tip sequence as the direct 
#'                    ancestor of the tree
#'                                                 
#' @return   A \code{graph} object representing the input tree. 
#'           
#' @details
#' Convert from phylo to graph object. Uses the node.label vector to label internal nodes. Nodes 
#' may rotate but overall topology will remain constant.
#' 
#' @references
#' \enumerate{
#'   \item  Hoehn KB, Lunter G, Pybus OG - A Phylogenetic Codon Substitution Model for Antibody 
#'              Lineages. Genetics 2017 206(1):417-427
#'              https://doi.org/10.1534/genetics.116.196303 
#'  \item  Hoehn KB, Vander Heiden JA, Zhou JQ, Lunter G, Pybus OG, Kleinstein SHK - 
#'              Repertoire-wide phylogenetic models of B cell molecular evolution reveal 
#'              evolutionary signatures of aging and vaccination. bioRxiv 2019  
#'              https://doi.org/10.1101/558825 
#' }
#'
#' @examples
#' \dontrun{
#'    library(igraph)
#'    library(ape)
#' 
#'    #convert to phylo
#'    phylo = graphToPhylo(graph)
#'    
#'    #plot tree using ape
#'    plot(phylo,show.node.label=TRUE)
#'    
#'    #store as newick tree
#'    write.tree(phylo,file="tree.newick")
#'    
#'    #read in tree from newick file
#'    phylo_r = read.tree("tree.newick")
#'    
#'    #convert to igraph
#'    graph_r = phyloToGraph(phylo_r,germline="Germline")
#'    
#'    #plot graph - same as before, possibly rotated
#'    plot(graph_r,layout=layout_as_tree)
#' }
#' 
#' @export
phyloToGraph <- function(phylo, germline=NULL) {
    names <- 1:length(unique(c(phylo$edge[, 1],phylo$edge[, 2])))
    for(i in 1:length(phylo$tip.label)){
        names[i] <- phylo$tip.label[i]
    }
    if(!is.null(phylo$node.label)){
        for(j in 1:length(phylo$node.label)){
            i <- i + 1
            names[i] <- phylo$node.label[j]
        }
    }
    d <- data.frame(cbind(phylo$edge,phylo$edge.length))
    names(d)=c("from", "to", "weight")

    if(!is.null(germline)){
        germnode <- which(phylo$tip.label == germline)
        phylo$uca = phylo$edge[phylo$edge[,2] == germnode,1]
        if(sum(d$from == phylo$uca) == 2){
            d[d$from == phylo$uca, ]$from <- germnode
            d <- d[!(d$from == germnode & d$to == germnode),] 
        }else{
            row <- which(d$from == phylo$uca & d$to == germnode)
            d[row,]$to <- phylo$uca
            d[row,]$from <- germnode
        }
    }

    d$to <- as.character(d$to)
    d$from <- as.character(d$from)
    g <- igraph::graph_from_data_frame(d)
    igraph::V(g)$name <- names[as.numeric(igraph::V(g)$name)]
    igraph::E(g)$label <- igraph::E(g)$weight
    return(g)
}


#' Convert a tree in igraph \code{graph} format to ape \code{phylo} format.
#' 
#' \code{graphToPhylo} a tree in igraph \code{graph} format to ape \code{phylo} 
#' format.
#' 
#' @param   graph  An igraph \code{graph} object.
#'
#' @return   A \code{phylo} object representing the input tree. Tip and internal node names are 
#'           stored in the \code{tip.label} and \code{node.label} vectors, respectively.
#'           
#' @details
#' Convert from igraph \code{graph} object to ape \code{phylo} object. If \code{graph} object
#' was previously rooted with the germline as the direct ancestor, this will re-attach the 
#' germline as a descendant node with a zero branch length to a new universal common ancestor (UCA) 
#' node and store the germline node ID in the \code{germid} attribute and UCA node number in 
#' the \code{uca} attribute. Otherwise these attributes will not be specified in the \code{phylo} object. 
#' Using \code{phyloToGraph(phylo, germline=phylo$germid)} creates a \code{graph} object with the germline 
#' back as the direct ancestor. Tip and internal node names are 
#' stored in the \code{tip.label} and \code{node.label} vectors, respectively.
#' 
#' @references
#' \enumerate{
#'   \item  Hoehn KB, Lunter G, Pybus OG - A Phylogenetic Codon Substitution Model for Antibody 
#'              Lineages. Genetics 2017 206(1):417-427
#'              https://doi.org/10.1534/genetics.116.196303 
#'  \item  Hoehn KB, Vander Heiden JA, Zhou JQ, Lunter G, Pybus OG, Kleinstein SHK - 
#'              Repertoire-wide phylogenetic models of B cell molecular evolution reveal 
#'              evolutionary signatures of aging and vaccination. bioRxiv 2019  
#'              https://doi.org/10.1101/558825 
#' }
#'
#' @examples
#' \dontrun{
#'    library(igraph)
#'    library(ape)
#' 
#'    #convert to phylo
#'    phylo = graphToPhylo(graph)
#'    
#'    #plot tree using ape
#'    plot(phylo,show.node.label=TRUE)
#'    
#'    #store as newick tree
#'    write.tree(phylo,file="tree.newick")
#'    
#'    #read in tree from newick file
#'    phylo_r = read.tree("tree.newick")
#'    
#'    #convert to igraph
#'    graph_r = phyloToGraph(phylo_r,germline="Germline")
#'    
#'    #plot graph - same as before, possibly rotated
#'    plot(graph_r,layout=layout_as_tree)
#' }
#' 
#' @export
graphToPhylo <- function(graph) {
    df  <- igraph::as_data_frame(graph)
    node_counts <- table(c(df$to,df$from))
    tips <- names(node_counts)[node_counts == 1]
    nodes <- names(node_counts)[node_counts > 1]

    germline <- tips[tips %in% df$from]
    if(length(germline) > 0){
        ucanode <- paste0(germline,"_UCA")#max(as.numeric(nodes))+1
        nodes <- c(ucanode,nodes)
        df[df$from == germline,]$from <- ucanode
        row <- c(ucanode,germline,0.0)
        names(row) <- c("from","to","weight")
        df <- rbind(df, row)
    }
    tipn <- 1:length(tips)
    names(tipn) <- tips
    noden <- (length(tips)+1):(length(tips)+length(nodes))
    names(noden) <- nodes
    renumber <- c(tipn,noden)

    df$from <- as.numeric(renumber[df$from])
    df$to <- as.numeric(renumber[df$to])    
    
    phylo <- list()
    phylo$edge <- matrix(cbind(df$from,df$to),ncol=2)
    phylo$edge.length <- as.numeric(df$weight)
    phylo$tip.label <- tips
    phylo$Nnode <- length(nodes)
    phylo$node.label <- nodes
    class(phylo) <- "phylo"

    if(length(germline) > 0){
        phylo <- rerootGermline(phylo, germline)
    }

    phylo = ape::ladderize(phylo, right=FALSE)
    
    return(phylo)
}

# Reroot phylogenetic tree to have its germline sequence at a zero-length branch 
# to a node which is the direct ancestor of the tree's UCA. Assigns \code{uca}
# to be the ancestral node to the tree's germline sequence, as \code{germid} as
# the tree's germline sequence ID. 
#
# @param   tree     An ape \code{phylo} object
# @param   germid   ID of the tree's predicted germline sequence
# @param   resolve  If \code{TRUE} reroots tree to specified germline sequnece.
#                   usually not necessary with IgPhyML trees analyzed with HLP model.
rerootGermline <- function(tree, germid, resolve=FALSE){
    if(resolve) {
        tree <- ape::root(phy=tree, outgroup=germid, resolve.root=T, edge.label=TRUE)
    }
    tree <- ape::reorder.phylo(tree, "postorder")  
    edges <- tree$edge
    rootnode <- which(tree$tip.label==germid)
    rootedge <- which(edges[, 2] == rootnode)
    rootanc <- edges[edges[, 2] == rootnode, 1]
    mrcaedge <- which(edges[, 1] == rootanc & edges[, 2] != rootnode)
    if(length(mrcaedge) > 1){
            print("POLYTOMY AT ROOT!")
            quit(save="no", status=1, runLast=FALSE)
    }
    tree$edge.length[mrcaedge] <- tree$edge.length[mrcaedge] + tree$edge.length[rootedge]
    tree$edge.length[rootedge] <- 0
    tree$uca <- rootanc
    tree$germid <- germid
    
    return(tree)
}

#' Read in output from IgPhyML
#' 
#' \code{readIgphyml} reads output from the IgPhyML phylogenetics inference package for 
#' B cell repertoires
#' 
#' @param    file          IgPhyML output file (.tab).
#' @param    id            ID to assign to output object.
#' @param    format        if \code{"graph"} return trees as igraph \code{graph} objects. 
#'                         if \code{"phylo"} return trees as ape \code{phylo} objects.
#' @param    branches        if \code{"distance"} branch lengths are in expected mutations per
#'                         site. If \code{"mutations"} branches are in expected mutations.
#' @param    collapse      if \code{TRUE} transform branch lengths to units of substitutions, 
#'                         rather than substitutions per site, and collapse internal nodes
#'                         separated by branches < 0.1 substitutions. Will also remove all
#'                         internal node labels, as it makes them inconsistent.
#'                                                
#' @return   A list containing IgPhyML model parameters and estimated lineage trees. 
#'           
#'           Object attributes:
#'           \itemize{
#'             \item  \code{param}:     Data.frame of parameter estimates for each clonal 
#'                                      lineage. Columns include: \code{CLONE}, which is the 
#'                                      clone id; \code{NSEQ}, the total number of sequences in 
#'                                      the lineage; \code{NSITE}, the number of codon sites;
#'                                      \code{TREE_LENGTH}, the sum of all branch lengths in 
#'                                      the estimated lineage tree; and \code{LHOOD}, the log 
#'                                      likelihood of the clone's sequences given the tree and
#'                                      parameters. Subsequent columns are parameter estimates 
#'                                      from IgPhyML, which will depend on the model used. 
#'                                      Parameter columns ending with \code{_MLE} are maximum 
#'                                      likelihood estimates; those ending with \code{_LCI} are 
#'                                      the lower 95%% confidence interval estimate; those ending 
#'                                      with \code{_UCI} are the upper 95%% confidence interval 
#'                                      estimate. The first line of \code{param} is for clone 
#'                                      \code{REPERTOIRE}, 
#'                                      which is a summary of all lineages within the repertoire.
#'                                      For this row, \code{NSEQ} is the total number of sequences, 
#'                                      \code{NSITE} is the average number of sites, and
#'                                      \code{TREE_LENGTH} is the mean tree length. For most 
#'                                      applications, parameter values will be the same for all 
#'                                      lineages within the repertoire, so access them simply by:
#'                                      \code{<object>$param$OMEGA_CDR_MLE[1]} to, for instance,
#'                                      get the estimate of dN/dS on the CDRs at the repertoire level.
#'             \item  \code{trees}:     List of tree objects estimated by IgPhyML. If 
#'                                      \code{format="graph"} these are igraph \code{graph} objects. 
#'                                      If \code{format="phylo"}, these are ape \code{phylo} objects.
#'             \item  \code{command}:   Command used to run IgPhyML.
#'           }
#'           
#' @details
#' \code{readIgphyml} reads output from the IgPhyML repertoire phylogenetics inference package. 
#' The resulting object is divded between parameter estimates (usually under the HLP19 model),
#' which provide information about mutation and selection pressure operating on the sequences.
#' 
#' Trees returned from this function are either igraph objects or phylo objects, and each may be 
#' visualized accordingly. Futher, branch lengths in tree may represent either the expected number of
#' substitutions per site (codon, if estimated under HLP or GY94 models), or the total number of 
#' expected substitutions per site. If the latter, internal nodes - but not tips - separated by branch
#' lengths less than 0.1 are collapsed to simplify viewing.
#' 
#' @references
#' \enumerate{
#'   \item  Hoehn KB, Lunter G, Pybus OG - A Phylogenetic Codon Substitution Model for Antibody 
#'              Lineages. Genetics 2017 206(1):417-427
#'              https://doi.org/10.1534/genetics.116.196303 
#'  \item  Hoehn KB, Vander Heiden JA, Zhou JQ, Lunter G, Pybus OG, Kleinstein SHK - 
#'              Repertoire-wide phylogenetic models of B cell molecular evolution reveal 
#'              evolutionary signatures of aging and vaccination. bioRxiv 2019  
#'              https://doi.org/10.1101/558825 
#' }
#'
#' @examples
#' \dontrun{
#'    # Read in and plot a tree from an igphyml run
#'    library(igraph)
#'    s1 <- readIgphyml("IB+7d_lineages_gy.tsv_igphyml_stats_hlp.tab", id="+7d")
#'    print(s1$param$OMEGA_CDR_MLE[1])
#'    plot(s1$trees[[1]], layout=layout_as_tree, edge.label=E(s1$trees[[1]])$weight)
#' }
#' 
#' @export
readIgphyml <- function(file, id=NULL, format=c("graph", "phylo"), collapse=FALSE,
    branches=c("distance","mutations")) {
    # Check arguments
    format <- match.arg(format)
    branches <- match.arg(branches)
    
    out <- list()
    trees <- list()
    df <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    params <- df[, !names(df) %in% c("TREE")]
    names(params) = tolower(names(params))
    out[["param"]] <- params
    out[["command"]] <- df[1, ]$TREE
    for (i in 2:nrow(df)) {
        tree <- ape::read.tree(text=df[i, ][["TREE"]])
        germ_base <- paste0(df[["CLONE"]][i], "_GERM")
        germ_id <- tree$tip.label[grepl(germ_base,tree$tip.label)]
        if(length(germ_id) > 1){
            stop("Can only be one tip of the form '<cloneid>_GERM'")
        }
        if(!ape::is.rooted(tree)){
            warning(paste("Tree",germ_base,"is not rooted and should be!"))
            tree <- rerootGermline(tree,germ_id,resolve=TRUE)
            tree$node.label <- NULL
        }
        rtree <- ape::ladderize(tree)
        rtree$germid <- germ_id
        if (branches == "mutations") {
            rtree$edge.length <- round(rtree$edge.length*df[i, ]$NSITE, digits=1)
        }
        if (collapse) {
            rtree$node.label <- NULL
            if(branches == "mutations"){
                rtree <- ape::di2multi(rtree, tol=0.1)
            }else{
                rtree <- ape::di2multi(rtree, tol=0.0001)
            }
        }
        if (format == "graph") {
            ig <- phyloToGraph(rtree, germline=rtree$germid)
            trees[[df[["CLONE"]][i]]] <- ig
        } else if (format == "phylo") {
            trees[[df[["CLONE"]][i]]] <- rtree
        } else {
            stop("Format must be either 'graph' or 'phylo'.")
        }
    }
    
    out[["trees"]] <- trees

    if (!is.null(id)) {
        out$param$id <- id
    }
    
    return(out)
}


#' Combine IgPhyML object parameters into a dataframe
#' 
#' \code{combineIgphyml} combines IgPhyML object parameters into a data.frame.
#' 
#' @param   iglist         list of igphyml objects returned by \link{readIgphyml}. 
#'                         Each must have an \code{id} column in its \code{param} attribute, 
#'                         which can be added automatically using the \code{id} option of 
#'                         \code{readIgphyml}.
#' @param   format         string specifying whether each column of the resulting data.frame
#'                         should represent a parameter (\code{wide}) or if 
#'                         there should only be three columns; i.e. id, varable, and value
#'                         (\code{long}).
#'                                                
#' @return   A data.frame containing HLP model parameter estimates for all igphyml objects.
#'           Only parameters shared among all objects will be returned.
#'           
#' @details
#' \code{combineIgphyml} combines repertoire-wide parameter estimates from mutliple igphyml
#' objects produced by readIgphyml into a dataframe that can be easily used for plotting and 
#' other hypothesis testing analyses.
#' 
#' All igphyml objects used must have an "id" column in their \code{param} attribute, which
#' can be added automatically from the \code{id} flag of \code{readIgphyml}. 
#' 
#' @references
#' \enumerate{
#'   \item  Hoehn KB, Lunter G, Pybus OG - A Phylogenetic Codon Substitution Model for Antibody 
#'              Lineages. Genetics 2017 206(1):417-427
#'              https://doi.org/10.1534/genetics.116.196303 
#'  \item  Hoehn KB, Vander Heiden JA, Zhou JQ, Lunter G, Pybus OG, Kleinstein SHK - 
#'              Repertoire-wide phylogenetic models of B cell molecular evolution reveal 
#'              evolutionary signatures of aging and vaccination. bioRxiv 2019  
#'              https://doi.org/10.1101/558825 
#' }
#'
#' @seealso  \link{readIgphyml} 
#'           
#' @examples
#' \dontrun{
#'    # Read in and combine two igphyml runs
#'    s1 <- readIgphyml("IB+7d_lineages_gy.tsv_igphyml_stats_hlp.tab", id="+7d")
#'    s2 <- readIgphyml("IB+7d_lineages_gy.tsv_igphyml_stats_hlp.tab", id="s2")
#'    combineIgphyml(list(s1, s2))
#' }
#' 
#' @export
combineIgphyml <- function(iglist, format=c("wide", "long")) {
    # Check arguments
    format <- match.arg(format)
    
    ordered_params <- c(
        "id", "nseq", "nsite", "lhood", "tree_length", 
        "omega_fwr_mle", "omega_fwr_lci", "omega_fwr_uci", 
        "omega_cdr_mle", "omega_cdr_lci", "omega_cdr_uci", 
        "kappa_mle", "kappa_lci", "kappa_uci", 
        "wrc_2_mle", "wrc_2_lci", "wrc_2_uci", 
        "gyw_0_mle", "gyw_0_lci", "gyw_0_uci", 
        "wa_1_mle", "wa_1_lci", "wa_1_uci", 
        "tw_0_mle", "tw_0_lci", "tw_0_uci", 
        "syc_2_mle", "syc_2_lci", "syc_2_uci", 
        "grs_0_mle", "grs_0_lci", "grs_0_uci")
    paramCount <- table(unlist(lapply(iglist, function(x) names(x$param))))
    params <- names(paramCount[paramCount == max(paramCount)])
    params <- ordered_params[ordered_params %in% params]
    if (sum(params == "id") == 0) {
        message <- "id not specified in objects. Use 'id' flag in readIgphyml."
        stop(message)
    }
    
    repertoires <- lapply(iglist, function(x) x$param[1, params])
    combined <- dplyr::bind_rows(repertoires)
    if (format == "long") {
        combined <- tidyr::gather(combined, "variable", "value", -!!rlang::sym("id"))
        combined$variable <- factor(combined$variable, levels=params)
    }
    
    return(combined)
}