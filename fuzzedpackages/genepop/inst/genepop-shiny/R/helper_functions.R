disk.usage <- function(path = Sys.getenv("HOME")) {
  ## FR-> system("which df", intern = TRUE) fails under windows because which is non-existent
  trywhich <- try(system("which df", intern = TRUE),silent=TRUE)
  if (inherits(trywhich,"try-error")) {
    message("Available disk space not checked (Windows?)")
    structure(c(Inf,Inf), names = c("used", "available"))
  } else {
    if(length(trywhich)) {
      cmd <- sprintf("df  %s", path)
      exec <- system(cmd, intern = TRUE)
      exec <- strsplit(exec[length(exec)], "[ ]+")[[1]]
      exec <- as.numeric(exec[3:4])
      structure(exec, names = c("used", "available"))
    } else {
      stop("'df' command not found")
    }
  }
}

genlab <-function (base, n)
{
    f1 <- function(cha, n) {
        if (nchar(cha) < n) {
            cha <- paste("0", cha, sep = "")
            return(f1(cha, n))
        }
        else {
            return(cha)
        }
    }
    w <- as.character(1:n)
    max0 <- max(nchar(w))
    w <- sapply(w, function(cha) f1(cha, max0))
    return(paste(base, w, sep = ""))
}

rmspaces <-function (charvec)
{
    charvec <- gsub("^([[:blank:]]*)([[:space:]]*)", "", charvec)
    charvec <- gsub("([[:blank:]]*)([[:space:]]*)$", "", charvec)
    return(charvec)
}

read.genepop <-function (file, ncode = 2L, quiet = FALSE) { # corrected from adegenet
    prevcall <- match.call()
    txt <- scan(file, sep = "\n", what = "character", quiet = TRUE)
    if (!quiet)
        cat("\nFile description: ", txt[1], "\n")
    txt <- txt[-1]
    txt <- gsub("\t", " ", txt)
    NA.char <- paste(rep("0", ncode), collapse = "")
    locinfo.idx <- 1:(min(grep("POP[ ]*$", toupper(txt))) - 1) ## line indices before the first pop
    locinfo <- txt[locinfo.idx] ## lines with locus names
    locinfo <- paste(locinfo, collapse = ",")
    loc.names <- unlist(strsplit(locinfo, "([,]|[\n])+"))
    loc.names <- rmspaces(loc.names)
    nloc <- length(loc.names)
    txt <- txt[-locinfo.idx]
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
    npop <- length(pop.idx)
    nocomma <- which(!(1:length(txt)) %in% grep(",", txt))
    splited <- nocomma[which(!nocomma %in% pop.idx)]
    if (length(splited) > 0) {
        for (i in sort(splited, decreasing = TRUE)) {
            txt[i - 1] <- paste(txt[i - 1], txt[i], sep = " ")
        }
        txt <- txt[-splited]
    }
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
    txt[length(txt) + 1] <- "POP"
    nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$",
        toupper(txt))) - 1
    pop <- factor(rep(1:npop, nind.bypop))
    txt <- txt[-c(pop.idx, length(txt))]
    temp <- sapply(1:length(txt), function(i) strsplit(txt[i],
        ","))
    ind.names <- sapply(temp, function(e) e[1])
    ind.names <- rmspaces(ind.names)
    vec.genot <- sapply(temp, function(e) e[2])
    vec.genot <- rmspaces(vec.genot)
    X <- matrix(unlist(strsplit(vec.genot, "[[:space:]]+")),
        ncol = nloc, byrow = TRUE)
    if (any(duplicated(ind.names))) {
        rownames(X) <- genlab("", nrow(X))
        warning("Duplicate individual names detected. Coercing them to be unique.")
    }
    else {
        rownames(X) <- ind.names
    }
    colnames(X) <- loc.names
    pop.names.idx <- cumsum(table(pop))
    pop.names <- ind.names[pop.names.idx]
    levels(pop) <- pop.names
    #if (!all(unique(nchar(X)) == (ncode * 2)))
    #    stop(paste("some alleles are not encoded with", ncode,
    #        "characters\nCheck 'ncode' argument"))
    #res <- df2genind(X = X, pop = pop, ploidy = 2, ncode = ncode,
    #    NA.char = NA.char)
    #res@call <- prevcall
    if (!quiet)
        cat("\n...done.\n\n")
    return(list(datas=X, pops=pop.names, nind.bypop=nind.bypop))
}
#read.genepop("~/projets/genepop/sample.txt")



get_html_content <- function (ficname)
{
    contenu = readLines(ficname)
    contenu = paste(contenu, sep="", collapse="</br>")
    return (contenu)
}

get_plots <- function(data)
{
  plots = vector("list", ncol(data)-1)
  st = "Status"
  for ( j in 1:(ncol(data)-1) )
  {
      colo  = colnames(data)[j]
      estet = aes_string(x = "Test", y = colo, fill = st )
      estet$x = ""
      plots[[j]] = ggplot(data, estet) +
      geom_bar(width = 1, stat = "identity", colour = "black") +
      coord_polar("y") +# scale_fill_brewer(palette="Dark2") +
      geom_text(aes_string(y = paste(colo,"/3 + c(0, cumsum(",colo,")[-length(",colo,")])", sep="" ), label=paste("percent(",colo,"/(sum(",colo,")))",sep="") ), position = position_jitter(h =0.1, w = 0.5))+ scale_fill_manual(values = c("#3366cc","white", "#ff9900" )) +
      labs(title=colo, x="", y="#cores")
  }
  return(plots)
}

getSeed <- function(seed){
  if(is.na(seed)) {
    return (12345678)
  } else {
    return (seed)
  }

}

getMSeed <- function(seed){
  if(is.na(seed)) {
    return (87654321)
  } else {
    return (seed)
  }

}

Gprint <- function(debug, message){
  if(debug){
    cat(message)
  }
}

DisableButtonGenepop <- function() {

  for (i in 1:8)  {
    for (j in 1:6) {
      boutton = paste("downloadOpt",i,j,"All",sep="")
      shinyjs::disable(boutton)}
  }
}

getFile <- function(file, dir="doc") {
  system.file(dir, file, package="genepop")
}
