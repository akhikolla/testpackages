# Do a regression to determine the optimal numbers of supported data to use; the minimum is 3
.check.equi <- function(dets, suggest=TRUE) {
  rawdata = dets[,4] #210Pb
  rawsd   = dets[,5] #sd(210pb)
  deps    = dets[,2] #depth

  lendat = length(rawdata)
  numdat = as.integer(.5*length(rawdata))
  usedat = rawdata[(lendat-3):lendat]
  usesd  = rawsd[(lendat-3):lendat]
  usex   = 1:length((lendat-3):lendat)
  usereg = lm(usedat ~ usex, weights=1/(usesd^2))
  reg    = coef(summary(usereg))[2,4]
  est    = coef(summary(usereg))[1,1]
  coe    = 3
  for(i in 1:numdat) {
    usedat = rawdata[(lendat-3-i):lendat]
    usesd  = rawsd[(lendat-3-i):lendat]
    usex   = 1:length((lendat-3-i):lendat)
    usereg = lm(usedat ~ as.numeric(scale(usex)), weights=1/(usesd^2))
    reg1   = coef(summary(usereg))[2,4]
    est1   = mean(usedat) #coef(summary(usereg))[1,1]
    if(reg1 > reg) {
      reg  = reg1
      coe  = (3+i)
      est  = est1
    }
  }
  
  if(suggest) {
    ans <- readline(message("The regression process proposes using the last ", as.integer(coe), " data points as estimates of the supported activity, with a p-value of ", round((reg),3), ", OK? (Y/n) "))
    if(!(ans=="y" || ans=="")) 
      stop("  OK. Please adapt settings.\n\n", call.=FALSE)
  } 
  c(coe, reg1)
}



# read the 210Pb dets file
.read.dets.plum <- function(core, coredir, n.supp=c(), date.sample, set=get('info'), sep=",", dec=".", cc=1, Bqkg = TRUE, radon.case=c(), suggest=TRUE) {

  # read the file
  csv.file <- paste0(coredir,  core, "/", core, ".csv")
  dat.file <- paste0(coredir,  core, "/", core, ".dat")

  changed <- 0 # if the file needs changes, it has to be written to at the end
  suggested.names <- c("labID","depth(cm)","density(g/cm^3)","210Pb(Bq/kg)","sd(210Pb)","thickness(cm)", "226Ra(Bq/kg)", "sd(226Ra)")
  if(file.exists(csv.file)) {
    dets <- read.table(csv.file, header=TRUE, sep=sep)
    if(file.exists(dat.file)) # deal with old .dat files
      if(file.info(csv.file)$mtime < file.info(dat.file)$mtime)
        message("Warning, the .dat file is newer than the .csv file! I will read the .csv file. From now on please modify ", csv.file, ", not ", dat.file, " \n", sep="") else
    message("Reading", csv.file, "\n")
  } else {
    if(file.exists(paste0(csv.file, ".txt"))) {
      file.rename(paste0(csv.file, ".txt"), csv.file)
      message("Removing .txt extension from .csv file")
    } else {
      message("No .csv file found, reading ", dat.file, " and converting it to .csv")
      dets <- read.table(dat.file, header=TRUE)
      changed <- 1
    }
  }

#    name <- tolower(names(dets))
  commas <- grep(",,", readLines(csv.file)) # check if there are too many commas (e.g., lines with just commas)
  if(length(!is.na(commas)) > 0) # often an artefact of spreadsheet programs
    stop("check the .csv file in a plain-text editor for 'orphan' commas.\n", call.=FALSE)

  # relations between the names of columns and their positions in the .csv file
  idColumn       = 1
  depthColumn    = 2
  rhoColumn      = 3 # density
  plumdataColumn = 4 # means of measurements
  stdColumn      = 5 # their errors
  deltaColumn    = 6 # sample thickess
  radonColumn    = 7 # if present
  sdRadonColumn  = 8 # if present

  date.infile = c(); nsupp.infile = c(); radoncase.infile = c()
  date.asoption = c(); nsupp.asoption = c(); radoncase.asoption = c()
  if(ncol(dets) == 6 || ncol(dets) == 8) # no additional information in file
    detsOrig = dets else
      if(ncol(dets) == 7 || ncol(dets) == 9) { # additional information in file
        n = ifelse(ncol(dets) == 7, 7, 9)
        detsOrig = dets[,-n]
      if(length(dets[1,n]) > 0)
        if(!is.na(dets[1,n]))
          date.infile = dets[1,n]
      if(length(dets[2,n]) > 0)
        if(!is.na(dets[2,n]))
          nsupp.infile = dets[2,n]
      if(length(dets[3,n]) > 0)
        if(!is.na(dets[3,n]))
          radoncase.infile = dets[3,n]
      } else
        stop(paste(csv.file, "should have between 6 and 9 columns. Please check."), call.=TRUE)

  # read sampling date
  if(length(date.infile) > 0)
    if(!is.numeric(date.infile)) # sampling date provided within the .csv file
      stop("The date (first number of last column of your .csv file) must be a numeric value.", call.=FALSE)
  if(length(date.sample) > 0) # then it's provided as option within the Plum command
    if(is.numeric(date.sample))
      date.asoption = date.sample else 
        stop("date.sample must be a numeric value.", call.=FALSE)

  # read n.supp (amount of tail measurements to estimate supported Pb210)
  nsupp.asoption = NA
  if(length(nsupp.infile) > 0)
    if(!is.numeric(nsupp.infile)) # sampling date provided within the .csv file
      stop("n.supp (second number of last column of your .csv file) must be a numeric value.", call.=FALSE)
  if(length(n.supp) > 0) # then provided as option within the Plum command
    if(is.numeric(n.supp))
      nsupp.asoption = n.supp else 
        stop("n.supp must be a numeric value.", call.=FALSE) 

  # read radon case
  radoncase.asoption = NA
  if(length(radoncase.infile) > 0)
    if(!is.numeric(radoncase.infile)) # sampling date provided within the .csv file
      stop("Radon.case (third number of last column of your .csv file) must be a numeric value.", call.=FALSE)
  if(length(radon.case) > 0) # then it's provided as option within the Plum command
    if(is.numeric(radon.case))
      radoncase.asoption = radon.case else 
        stop("radon.case must be a numeric value.", call.=FALSE) 
  
  # determine sampling date
  if(length(date.infile) == 0) {
    if(length(date.asoption) == 0) {
      ans = readline("Please provide a date (in AD) for when the Pb210 samples were measured: ")
      if( grepl("^[0-9.][0-9]*[.]?[0-9]*[0-9.]$",ans) == FALSE ) 
        if( grepl("^[0-9]+$", ans) == FALSE )
        #  if(!is.numeric(ans))  
          stop("date.sample must be a numeric value, e.g., 2019.5.", call.=FALSE)
      date.sample = as.numeric(ans)
    } else
      date.sample = date.asoption
  } else {
     if(length(date.asoption) == 0) {
       message("Using date.sample provided in .csv file, ", date.infile, ".")
       date.sample = date.infile
     } else {
       message("date.sample provided both in the .csv file and as option; using the one from the file, ", date.infile, ".")
       date.sample = date.infile
    }
  }

  # determine n.supp. Only needed if no radon provided - could be used if radon provided but assuming constant supported Pb.
  if(length(nsupp.infile) == 0) { 
    if(length(nsupp.asoption) == 0 || is.na(nsupp.asoption)) 
      n.supp = .check.equi( detsOrig, suggest=suggest )[1]  else 
        n.supp = nsupp.asoption
    } else {
        n.supp = nsupp.infile # infile declaration of n.supp takes precedence over n.supp as option
        if(length(nsupp.asoption) > 0)
          message("Both the .csv file and the n.supp option provide n.supp; using the one from the file, ", nsupp.infile, ".")
    }

# determine radon.case
  if(length(radoncase.infile) == 0) {
    if(length(radoncase.asoption) == 0 || is.na(radoncase.asoption)) {
      if(ncol(detsOrig) == 6) {
        message("No radon present, radon case not given, setting it at 0.")
        radon.case = 0
      } else 
           stop("Please specify radon.case=1 or radon.case=2 (in .csv file or as Plum option)", call.=FALSE)
    } else 
         radon.case = radoncase.asoption
    } else {
        radon.case = radoncase.infile
        if(length(radoncase.asoption) == 0)
          message("Using radon case ", radoncase.infile, ifelse(radoncase.infile < 2, " (constant", " (varying"), " supported Pb).") else
            message("Both the .csv file and the radon.case option provide radon.case; using the one from the file, ", radoncase.infile, ".")
      }    

# now check for e.g., absence of radon, ...
  if(radon.case == 0) {
    if(ncol(detsOrig) > 6) {
      message("Found radon.case = 0 but there is radon. Therefore, setting radon.case to 1 (constant supported Pb).")
      radon.case = 1
      }
  }
  if(radon.case == 1) 
    if(ncol(detsOrig) == 6) {
      message("Found radon.case = 1 but there is no radon. Therefore, setting radon.case to 0 (constant supported Pb).")
      radon.case = 0
  }
  if(radon.case == 2) {
    if(ncol(detsOrig) == 6)
       stop("Cannot have radon.case 2 (varying supported Pb-210) without radon measurements.", call. = FALSE)
    if(ncol(detsOrig) == 8) {
      if(n.supp > 0) {
        message("The radon case can't be 2 when the number of supported measurements is >0; setting n.supp to 0.")
        n.supp = 0
      }
    }
  }
 
  #check that depths are in ascending order
  if(min(diff(dets[,depthColumn])) < 0) {
    message("Warning, the depths are not in ascending order, I will correct this.")
    dets <- dets[ order(dets[,depthColumn]),]
    write.table(dets, csv.file, sep=sep, dec=dec, row.names=FALSE, quote=FALSE)
  }


  if( core == "HP1C" ) {
    radonColumn = 4
    sdRadonColumn = 5
    supportedData = detsOrig[30:33,c(radonColumn, sdRadonColumn)] # so, removed the bottommost few samples
    detsOrig = detsOrig[1:29,]
    dets = dets[1:29,]
    if( radon.case < 0 )
      radon.case = 0
    date.sample = 2018.5
  } else 
     if( ncol(detsOrig) == 6 ) { # then repeat the Pb210 columns again
       radonColumn = 4
       sdRadonColumn = 5

    if(length(n.supp) == 0 || is.na(n.supp) || n.supp == 0 ) { # the number of supported data
      if( radon.case < 0 )
        radon.case = 0 # no radon, estimate supported data from the tail where 210Pb reaches equilibrium
      supportedData <- c()
    } else if( n.supp > 0 ) {
      if( radon.case < 0 )
        radon.case = 1 # assuming constant supported radon
      supportedData = detsOrig[(nrow(detsOrig)-n.supp+1):nrow(detsOrig),c(radonColumn, sdRadonColumn)]
      # n.supp <<- n.supp
      detsOrig = detsOrig[1:(nrow(detsOrig)-n.supp),]
    } else { # estimate supported Pb210 from the tail measurements, using linear regression
      tmp = .check.equi(detsOrig)
      n.supp = tmp[1]

      supportedData = detsOrig[(nrow(detsOrig)-n.supp+1):nrow(detsOrig),c(radonColumn, sdRadonColumn)]
  #    dets = dets[1:(nrow(dets)-n.supp),]
      if( radon.case < 0 )
        radon.case = 1
    }
  } else if( ncol(detsOrig)  == 8 ) {
    radonColumn    = 7
    sdRadonColumn  = 8
    #columns 7 and 8 are supported data
    supportedData = detsOrig[,c(7, 8)]
    detsOrig = detsOrig[,-c(radonColumn,sdRadonColumn)]

    if( length(is.na(supportedData)) > 0 ) {
      message("Missing values are detected; the radon case is set to 1.")

      elim <- c() # get rid of data with NAs
      for(i in 1:nrow(supportedData)) 
        if( length(is.na(supportedData[i,])) > 0 ) 
          elim <- c( elim, i )
      supportedData = supportedData[ -elim, ]
      
    #  supportedData <<- supportedData

      radon.case = 1
      ans <- readline(message("Additionally, do you want to set a number of tail measurements for supported Pb-210? (y/n)"))
      if(tolower(substr(ans, 1, 1)) == "y") {
        ans <- readline(message("Ok, n.supp="))
        n.supp = as.integer(ans)
        radonColumn = 4
        sdRadonColumn = 5
        tmp <- detsOrig[(nrow(detsOrig)-n.supp+1):(nrow(detsOrig)),c(radonColumn, sdRadonColumn)]
        names(tmp) <- colnames(supportedData)
        supportedData <- rbind( supportedData,  tmp)
        detsOrig <- detsOrig[1:(nrow(detsOrig)-n.supp),]
      }

    } else if( n.supp > 0 ) {
      if( radon.case < 0 )
        radon.case = 1
      radonColumn = 4
      sdRadonColumn = 5
      tmp <- detsOrig[(nrow(detsOrig)-n.supp+1):(nrow(detsOrig)),c(radonColumn, sdRadonColumn)]
      names(tmp) <- colnames(supportedData)
      supportedData <- rbind( supportedData,  tmp)
      detsOrig <- detsOrig[1:(nrow(detsOrig)-n.supp),]
    } else 
       if(suggest) { # was !provided
          message("Plum can assume to have a constant supported 210Pb and use the 226Ra data to infer this one value.")
          message(" Alternatively, we can also assume individual supported 210Pb values per measured depth.")
          message(" It is important to consider that this will greatly increase the computing time and it should only be used when clear patterns are observed in the 226Ra data.")
          ans <- readline(message(" Do you want to use the individual supported 210Pb? (y/n) "))
          if( !tolower(substr(ans, 1, 1)) == "y") {
            message(" OK, assuming constant supported 210Pb\n")
           radon.case = 1
         } else {
           message(" OK, using individual supported 210Pb per data point.")
           radon.case = 2
        }
    }
  } else {
    stop("Unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct dets file.", call.=FALSE)
  }

  # more sanity checks
  if(!is.numeric(dets[,plumdataColumn]) || !is.numeric(dets[,stdColumn]) || !is.numeric(dets[,depthColumn]))
    stop("unexpected values in dets file, I expected numbers. Check the manual.", call.=FALSE)
  if(!is.numeric(dets[,deltaColumn]) || !is.numeric(dets[,rhoColumn]) )
    stop("unexpected values in dets file, I expected numbers. Check the manual.", call.=FALSE)


  # if current dets differ from original .csv file, rewrite it
  if(changed > 0)
    write.table(dets, csv.file, sep=paste(sep, "\t", sep=""), dec=dec, row.names=FALSE, col.names=suggested.names[1:ncol(dets)], quote=FALSE)

  dets = dets[,c(idColumn, plumdataColumn, stdColumn, depthColumn, deltaColumn, rhoColumn)]

  # find the plot limits
  if(ncol(detsOrig) == 6) {
    age.min = min( c(detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]-detsOrig[,5]) )
    age.max = max( c(detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]+detsOrig[,5]) )
  } else {
    age.min = min( c(detsOrig[,2]-(detsOrig[,6]/2),detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]-detsOrig[,5],detsOrig[,7]-detsOrig[,8]) )
    age.max = max( c(detsOrig[,2]-(detsOrig[,6]/2),detsOrig[,2]-(detsOrig[,6]/2)), c( detsOrig[,4]+detsOrig[,5],detsOrig[,7]+detsOrig[,8]) )
  }

  # plot the data
  layout(1)
  oldpar <- par(mar=c(3,3,1,1), mgp=c(1.5,.7,.0), bty="l")
  on.exit(par(oldpar))

  age.lim = extendrange(c(age.min, age.max), f=0.01)
  dlim = c(0, max(detsOrig[,depthColumn]))
  ylab <- ifelse(Bqkg, '210Pb (Bq/kg)', '210Pb (dpm/g)')
  plot(0, type='n', pch=16,col=c(rep('red',nrow(detsOrig)),rep('red',nrow(detsOrig))),
    cex=.3, ylab=ylab, xlab='depth(cm)', xlim = dlim, ylim = age.lim )

 detsOrig <<- detsOrig
    
  rect(detsOrig[,2]-detsOrig[,6], detsOrig[,4]-detsOrig[,5],
    detsOrig[,2], detsOrig[,4]+detsOrig[,5],
    lty=3, border=4)
  if(ncol(detsOrig) > 6)
    rect(detsOrig[,2], detsOrig[,7]-detsOrig[,8], 
      detsOrig[,2]-detsOrig[,6], detsOrig[,7]+detsOrig[,8],
      lty=3, border=2)
      

  list(dets, supportedData, radon.case, date.sample, detsOrig, n.supp)
}



#' @name Plum.cleanup
#' @title Remove files made to produce the current core's age-depth model.
#' @description Remove files .bacon, .out, .pdf, _ages.txt, and _settings.txt of current core.
#' @details If cores behave badly, you can try cleaning up previous runs and settings, by
#' removing files .bacon, .out, .pdf, _ages.txt, and _settings.txt of current core.
#' @return A message stating that the files and settings of this run have been deleted.
#' @param set Detailed information of the current run, stored within this session's memory as variable \code{info}.
#' @author Maarten Blaauw, J. Andres Christen
#' @seealso \url{http://www.qub.ac.uk/chrono/blaauw/manualBacon_2.3.pdf}
#' @export
Plum.cleanup <- function(set=get('info')) {
  files <- c(paste0(set$prefix, ".bacon"), paste0(set$prefix, ".out"),
    paste0(set$prefix, ".pdf"), paste0(set$prefix, "_ages.txt"),
    paste0(set$coredir,set$core, "/", set$core, "_settings.txt"))
  for(i in files)
    if(file.exists(i))
      tmp <- file.remove(i)
  if(exists("tmp"))
    rm(tmp)
  message("Previous Plum runs of core", set$core, "with thick =", set$thick, "deleted. Now try running the core again\n")
}



# read in default values, values from previous run, any specified values, and report the desired one. Internal function.
.plum.settings <- function(core, coredir, dets, thick, remember=TRUE, d.min, d.max, d.by, depths.file,
  slump, acc.mean, acc.shape, mem.mean, mem.strength, boundary, hiatus.depths, hiatus.max, hiatus.shape,
  BCAD, cc, postbomb, cc1, cc2, cc3, cc4, depth.unit, normal, t.a, t.b, delta.R, delta.STD, prob,
  defaults, runname, ssize, dark, MinAge, MaxAge, cutoff, age.res, after, age.unit,
  supportedData, date.sample, Al, phi.shape, phi.mean, s.shape, s.mean, radon.case, Bqkg, n.supp) {

  vals <- list(d.min, d.max, d.by, depths.file, slump, acc.mean, acc.shape, mem.mean, mem.strength, boundary, hiatus.depths, hiatus.max, BCAD, cc, postbomb, cc1, cc2, cc3, cc4, depth.unit, normal, t.a, t.b, delta.R, delta.STD, prob, age.unit)
  valnames <- c("d.min", "d.max", "d.by", "depths.file", "slump", "acc.mean", "acc.shape", "mem.mean", "mem.strength", "boundary", "hiatus.depths", "hiatus.max", "BCAD", "cc", "postbomb", "cc1", "cc2", "cc3", "cc4", "depth.unit", "normal", "t.a", "t.b", "delta.R", "delta.STD", "prob", "age.unit")
  #TODO: modificar para que acepte el vector de soportado y los valores propios de plum, como es el nombre de archivo de soportado y los parametros como "Al"
  extr <- function(i, def=deffile, pre=prevfile, exists.pre=prevf, rem=remember, sep=" ", isnum=TRUE) {
    if(length(vals[[i]]) > 0) # tmp
      if(any(is.na(vals[[i]]))) {
        ext.def <- strsplit(def[i], sep)[[1]]
        ext.def <- ext.def[-length(ext.def)] # remove description
        if(exists.pre) {
          ext.pre <- strsplit(pre[i], sep)[[1]]
          ext.pre <- ext.pre[-length(ext.pre)] # remove description
          if(def[i] == pre[i]) # values for dev and pre similar, no worries
            ext <- ext.pre else
              if(rem) {
                if(i==13) ifelse(ext.pre, "using BC/AD", "using cal BP") else
                if(i>2) message(" using previous run's value for ", valnames[i], ", ", ext.pre)
                ext <- ext.pre
              } else {
                  if(i==13) ifelse(ext.def, "using BC/AD", "using cal BP") else
                  if(i>2) message(" using default value for ", valnames[i], ", ", ext.def)
                  ext <- ext.def
                }
        } else ext <- ext.def

        if(any(ext=="NA") || any(is.na(ext))) NA else
          if(isnum) as.numeric(ext) else noquote(ext)
      } else
        if(isnum) as.numeric(vals[[i]]) else vals[[i]]
  }

  # read in default values and those of previous run if available
  deffile <- readLines(defaults, n=-1)
  prevfile <- paste(coredir, core, "/", core, "_settings.txt", sep="")
  prevf <- FALSE
  if(file.exists(prevfile)) {
    prevfile <- readLines(prevfile, n=-1)
    if(length(prevfile) > 0) prevf <- TRUE
  }

  #d.min <- extr(1); d.by <- extr(3); depths.file <- extr(4)
  #slump <- extr(5); acc.mean <- extr(6);
  #if(length(acc.shape) == 1)
  #  acc.shape <- extr(7)
  #mem.mean <- extr(8)
  #mem.strength <- extr(9)
  #boundary <- if(is.na(boundary)[1]) NA else sort(extr(10))
  #hiatus.depths <- if(is.na(hiatus.depths)[1]) NA else sort(extr(11))
  #hiatus.max <- extr(12)
  #BCAD <- extr(13); cc <- extr(14); postbomb <- extr(15); cc1 <- extr(16, isnum=FALSE)
  #cc2 <- extr(17, isnum=FALSE); cc3 <- extr(18, isnum=FALSE); cc4 <- extr(19, isnum=FALSE)
  #depth.unit <- extr(20, isnum=FALSE); normal <- extr(21); t.a <- extr(22); t.b <- extr(23)
  #delta.R <- extr(24); delta.STD <- extr(25); prob <- extr(26); age.unit <- extr(27, isnum=FALSE)

  if(is.na(d.min) || d.min=="NA")
    d.min <- min(dets[,4])
  if(is.na(d.max) || d.max=="NA")
    d.max <- max(dets[,4])
  if(length(acc.shape) < length(acc.mean))
    acc.shape <- rep(acc.shape, length(acc.mean)) else
      if(length(acc.shape) > length(acc.mean))
        acc.mean <- rep(acc.mean, length(acc.shape))
  if(length(mem.strength) < length(mem.mean))
    mem.strength <- rep(mem.strength, length(mem.mean)) else
      if(length(mem.strength) > length(mem.mean))
        mem.mean <- rep(mem.mean, length(mem.strength))

  ## produce/update settings file, and return the values
  prevfile <- file(paste(coredir, core, "/", core, "_settings.txt", sep=""), "w")
  scat <- function(m, n="") cat(m, n, sep="", file=prevfile)
  cat(d.min, " #d.min\n", d.max, " #d.max\n", d.by, " #d.by\n",
    depths.file, " #depths.file\n", slump, " #slump\n", sep="", file=prevfile)
  for(i in acc.mean) scat(i, " "); scat("#acc.mean\n")
  for(i in acc.shape) scat(i, " "); scat("#acc.shape\n", "")
  for(i in mem.mean) scat(i, " "); scat("#mem.mean\n", "")
  for(i in mem.strength) scat(i, " "); scat("#mem.strength\n", "")
  for(i in boundary) scat(i, " "); scat("#boundary\n", "")
  for(i in hiatus.depths) scat(i, " "); scat("#hiatus.depths\n", "")
  for(i in hiatus.max) scat(i, " "); scat("#hiatus.max\n", "")
  #for(i in hiatus.max) scat(i, " "); scat("#hiatus.max\n", "") # redundant
  cat(BCAD, " #BCAD\n", cc, " #cc\n", postbomb, " #postbomb\n",
    cc1, " #cc1\n", cc2, " #cc2\n", cc3, " #cc3\n", cc4, " #cc4\n",
    depth.unit, " #depth.unit\n", normal, " #normal\n", t.a, " #t.a\n", t.b, " #t.b\n",
    delta.R, " #delta.R\n", delta.STD, " #d.STD\n", prob, " #prob\n", age.unit, "#age.unit\n", sep="", file=prevfile)

  cat(date.sample, " #date.sample\n", Al, " #Al\n", phi.shape, " #phi.shape\n", phi.mean, " #phi.mean\n",
    s.shape, " #s.shape\n", s.mean, " #s.mean\n", radon.case, " #radon.case\n", Bqkg, " #Bqkg\n", sep="", file=prevfile)

  cat(n.supp, " #n.supp\n", sep="", file=prevfile);

  close(prevfile)

  if(length(MinAge) == 0)
    MinAge <- min(1950 - as.integer(format(Sys.time(), "%Y")), round(dets[,2] - (5*dets[,3])))
  if(length(MaxAge) == 0)
    MaxAge <- max(1e6, round(dets[,2] + (5*dets[,3])))

  theta0 = 1950 - date.sample

  list(core=core, thick=thick, dets=dets, d.min=d.min, d.max=d.max, coredir=core,
    d.by=d.by, depths.file=depths.file, slump=slump,
    acc.mean=acc.mean, acc.shape=acc.shape, mem.mean=mem.mean,
    mem.strength=mem.strength, boundary=boundary,
    hiatus.depths=hiatus.depths, hiatus.max=hiatus.max,
    BCAD=BCAD, cc=cc, postbomb=postbomb,
    cc1=cc1, cc2=cc2, cc3=cc3, cc4=cc4, depth.unit=noquote(depth.unit), unit=depth.unit, age.unit=noquote(age.unit), normal=normal,
    t.a=t.a, t.b=t.b, delta.R=delta.R, delta.STD=delta.STD, prob=prob, date=date(),
    runname=runname, ssize=ssize, dark=dark, MinAge=MinAge, MaxAge=MaxAge,
    cutoff=cutoff, age.res=age.res, after=after,
    supportedData=supportedData, theta0 = theta0, Al=Al, phi.shape=phi.shape, phi.mean=phi.mean, s.shape=s.shape, s.mean=s.mean,
    radon.case=radon.case, Bqkg=Bqkg)
}



#function to merge dets of plum and bacon data
.merge.dets <- function(detsPlum, detsBacon, delta.R, delta.STD, t.a, t.b, cc){
  if( ncol(detsBacon) >= 5 ){
    cc <- detsBacon[,5]
    detsBacon <- detsBacon[,-5]
  }else{
    cc <- array(cc, dim=c(nrow(detsBacon),1))
  }

  if( ncol(detsBacon) < 9 ){

    for(i in (ncol(detsBacon)+1):9){
      if( i == 5){
        col <- array(delta.R, dim=c(nrow(detsBacon),1))
      }else if(i == 6){
        col <- array(delta.STD, dim=c(nrow(detsBacon),1))
      }else if(i == 7){
        col <- array(t.a, dim=c(nrow(detsBacon),1))
      }else if(i == 8){
        col <- array(t.b, dim=c(nrow(detsBacon),1))
      }else if(i==9){
        col <- cc
      }
      detsBacon <- cbind(detsBacon, col)
    }
    colnames(detsBacon) <- c("labID", "X210Pb.Bq.kg.", "sd.210Pb.", "depth.cm.", "thickness.cm.", "density.g.cm.3.",  "t.a", "t.b", "cc")
    #print(detsBacon)
  }

  if( ncol(detsPlum) < 9 ){
    for(i in (ncol(detsPlum)+1):9){
      if( i == 5){
        col <- array(delta.R, dim=c(nrow(detsPlum),1))
      }else if(i == 6){
        col <- array(delta.STD, dim=c(nrow(detsPlum),1))
      }else if(i == 7){
        col <- array(t.a, dim=c(nrow(detsPlum),1))
      }else if(i == 8){
        col <- array(t.b, dim=c(nrow(detsPlum),1))
      }else if(i==9){
        col <- array(5, dim=c(nrow(detsPlum),1))
      }
      detsPlum <- cbind(detsPlum, col)
    }
    colnames(detsPlum) <- c("labID", "X210Pb.Bq.kg.", "sd.210Pb.", "depth.cm.", "thickness.cm.", "density.g.cm.3.",  "t.a", "t.b", "cc")
    #print(detsPlum)
  }

  dets <- rbind(detsPlum, detsBacon, make.row.names = FALSE)
  dets <- dets[ order(dets[,4]),]

}



# write files to be read by the main Bacon age-depth modelling function
.write.plum.file <- function(set=get('info')) {

  #a relation between the name of column and his position
  #These are the column of the plum file
  idColumn       = 1
  plumdataColumn = 2
  stdColumn      = 3
  depthColumn    = 4
  deltaColumn    = 5
  rhoColumn      = 6

  if(length(set$slump) > 0) {
    dets <- set$slumpdets
    hiatus.depths <- set$slumphiatus
    boundary <- set$slumpboundary
  } else {
    dets <- set$dets
    hiatus.depths <- set$hiatus.depths
    boundary <- set$boundary
  }

  if(is.na(set$d.min) || set$d.min < min(dets[,depthColumn])) { # repeat relevant row, change error and depth
    # extrap <- c(NA, min(dets[,2]), max(1e5, 1e3*dets[,2], 1e3*dets[,3]), set$d.min, 0)
    dets <- rbind(dets[which(dets[,depthColumn] == min(dets[,depthColumn]))[1],], dets, make.row.names=FALSE)
    dets[1,1] <- NA # calling this "d.min" causes issues
    dets[1,3] <- max(1e5, 1e3*dets[,4], 1e3*dets[,3])
    dets[1,depthColumn] <- set$d.min
  }
  #print(set$d.max)
  #print(max(dets[,depthColumn]))
  if(is.na(set$d.max) || set$d.max > max(dets[,depthColumn])) { # repeat relevant row, change error and depth
    # extrap <- c(NA, max(dets[,2]), max(1e5, 1e3*dets[,2], 1e3*dets[,3]), set$d.max, 0)
    dets <- rbind(dets, dets[which(dets[,depthColumn] == max(dets[,depthColumn]))[1],], make.row.names=FALSE)
    dets[nrow(dets),1] <- NA # calling this "d.max" causes issues
    dets[nrow(dets),3] <- max(1e5, 1e3*dets[,4], 1e3*dets[,3])
    dets[nrow(dets),depthColumn] <- set$d.max
  }

  #print(dets)

  supportedData <- set$supportedData

  fl <- file(set$bacon.file, "w")
  cat("## Ran on", set$date, "\n\n", file=fl)
  cat("Cal 0 : ConstCal;\nCal 1 : ",
  if(set$cc1=="IntCal20" || set$cc1=="\"IntCal20\"") "IntCal20"
    else noquote(set$cc1), ", ", set$postbomb, ";\nCal 2 : ",
  if(set$cc2=="Marine20" || set$cc2=="\"Marine20\"") "Marine20"
    else noquote(set$cc2), ";\nCal 3 : ",
  if(set$cc3=="SHCal20" || set$cc3=="\"SHCal20\"") "SHCal20"
    else noquote(set$cc3), ", ", set$postbomb, ";",
  if(set$cc4=="ConstCal" || set$cc4=="\"ConstCal\"") set$cc4 <- c()
    else
      paste("\nCal 4 : GenericCal, ", set$cc4, ";", sep=""), sep="", file=fl)
  cat("\nCal 4 : ConstCal;", sep="", file=fl)

  cat("\n##          alPhi mPhi  alS  mS     Al   theta0  Radon_case  supported_data_file", file=fl)

  cat("\nCal 5 : Plum, ", set$phi.shape, ", ",  set$phi.mean, ", ",  set$s.shape, ", ", set$s.mean, ", ", set$Al, ", ", set$theta0, ", ",
        set$radon.case, ", ", set$plum.file,";", sep="", file=fl)

  #cat("\n##   id.    210Pb   std   depth   delta     rho  t.a t.b cc ... Plum: 210Pb data", file=fl)
  cat("\n##    ", colnames(dets), " ... Plum: 210Pb data",sep=", ", file=fl)

  # we need to send the dets with all columns so pre-processing is needed
  for( i in 1:nrow(dets) ){
    cat( "\nDet ", i-1, " : ", as.character(dets[i,1]),
        " , ", dets[i,2],
        ", ", dets[i,3],
        ", ", dets[i,4],
        ", ", dets[i,5],
        ", ", dets[i,6],
        ", ", dets[i,7],
        ", ", dets[i,8],
        ", ", dets[i,9],
        ";", sep="", file=fl)
  }

  if(!is.na(hiatus.depths[1])) {
    if(is.null(boundary[1]))
      message("\n  Hiatus set at depth(s)", hiatus.depths, "\n") else
        message("\n  Boundary set at depth(s)", boundary, "\n")
    if(length(set$acc.shape)==1)
      set$acc.shape <- rep(set$acc.shape, length(hiatus.depths)+1)
    if(length(set$acc.mean)==1)
      set$acc.mean <- rep(set$acc.mean, length(hiatus.depths)+1)
    if(length(set$hiatus.max)==1)
      set$hiatus.max <- rep(set$hiatus.max, length(hiatus.depths))
#      if(length(set$hiatus.shape)==1)
#        set$hiatus.shape <- rep(set$hiatus.shape, length(set$hiatus.depths))
    .assign_to_global("info", set)
    cat("\n\n### Depths and priors for fixed hiatuses, in descending order",
      "\n##### cm  alpha beta      ha     hb", file=fl)
    for(i in length(hiatus.depths):1)
      cat("\nHiatus ", i-1, ":  ", hiatus.depths[i], ",  ", set$acc.shape[i+1],
        ",  ", set$acc.shape[i+1]/set$acc.mean[i+1], ",  ", .1, # last value (h.a) was NA but this conflicts with setting initial values for hiatus length
        ",  ", set$hiatus.max[i], ";", sep="", file=fl)
  }

  cK <- set$d.min+(set$thick*set$K)
  ### final parameters - dmax now calculated as dmin+(dC*K)
  if( is.na(set$seed) ){
  wrapup <- paste("\n\n##\t\t K   MinAge   MaxAge   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax",
    "\nBacon 0: ", ifelse(set$normal, "FixNor", "FixT"), ", ", set$K,
    ",  ", set$theta0-.02, ",  ", 26500, ",  ", set$theta0-0.01, ",  ", set$theta0+0.01,
    ",  ", set$mem.strength*set$mem.mean, ",  ", set$mem.strength*(1-set$mem.mean),
    ",  ", set$acc.shape[1], ",  ", set$acc.shape[1]/set$acc.mean[1], ", ", set$d.min,
    ", ", cK, ";\n", sep="")
  }else{
    wrapup <- paste("\n\n##\t\t K   MinAge   MaxAge   th0   th0p   w.a   w.b   alpha  beta  dmin  dmax  seed",
      "\nBacon 0: ", ifelse(set$normal, "FixNor", "FixT"), ", ", set$K,
      ",  ", set$theta0-.02, ",  ", 26500, ",  ", set$theta0-0.01, ",  ", set$theta0+0.01,
      ",  ", set$mem.strength*set$mem.mean, ",  ", set$mem.strength*(1-set$mem.mean),
      ",  ", set$acc.shape[1], ",  ", set$acc.shape[1]/set$acc.mean[1], ", ", set$d.min,
      ", ", cK, ", ", set$seed, ";\n", sep="")
  }
  cat(wrapup, file=fl)
  close(fl)

  fl <- file(set$plum.file, "w")
  for (i in 1:nrow(supportedData)){
    for (j in 1:ncol(supportedData)){
      cat( supportedData[i,j], " ", sep="", file= fl )
    }
    cat("\n", file=fl)
  }
  close(fl)
}



# function to read output files into memory
.Plum.AnaOut <- function(fnam, set=get('info')) {
  out <- read.table(fnam)
  n <- ncol(out)-1
  set$nPs  <- n
  set$TrPs <- nrow(out)
  set$phi  <- out[,1]
  set$ps   <- out[,2:(n+1)]
  set
}



# read the dets file, converting old formats to new ones if so required
.read.dets.plumbacon <- function(core, otherdates, coredir, set=get('info'), sep=",", dec=".", cc=1) {
  # if a .csv file exists, read it (checking that it is more recent than any .dat file in the folder). Otherwise, read the .dat file, check the columns, report back if >4 (>5?) columns, and convert to .csv (report this also)
  
  dat.file <- paste0(coredir,  core, "/", otherdates, ".dat")
  if(length(grep(".csv", otherdates)) > 0) # if the name has extension .csv 
  	csv.file = paste0(coredir, core, "/", otherdates) else
      csv.file <- paste0(coredir,  core, "/", otherdates, ".csv")

  dR.names <- c("r", "d", "d.r", "dr", "deltar", "r.mn", "rm", "rmn", "res.mean", "res.mn", "delta.r")
  dSTD.names <- c("d.std", "std", "std.1", "dstd", "r.std", "rstd", "res.sd", "delta.std", "deltastd")
  ta.names <- c("t", "t.a", "ta", "sta")
  tb.names <- c("t", "t.b", "tb", "stb")
  cc.names <- c("c", "cc")
  suggested.names <- c("labID", "age", "error", "depth", "cc", "dR", "dSTD", "ta", "tb")
  changed <- 0

  if(file.exists(csv.file)) {
    dets <- read.table(csv.file, header=TRUE, sep=sep)
    if(file.exists(dat.file)) # deal with old .dat files
      if(file.info(csv.file)$mtime < file.info(dat.file)$mtime)
        message("Warning, the .dat file is newer than the .csv file! I will read the .csv file. From now on please modify ", csv.file, ", not ", dat.file) else
          message("Reading", csv.file)
    } else {
      if(file.exists(paste0(csv.file, ".txt"))) {
        file.rename(paste0(csv.file, ".txt"), csv.file)
        message("Removing .txt extension from .csv file")
      } else {
        message("No .csv file found, reading", dat.file, "and converting it to .csv")
        dets <- read.table(dat.file, header=TRUE)
        changed <- 1
        }
    }
  name <- tolower(names(dets))
  commas <- grep(",,", readLines(csv.file)) # check if there are too many commas (e.g., lines with just commas)
  if(length(!is.na(commas)) > 0) # often an artifact of spreadsheet programs
    stop("check the .csv file in a plain-text editor for 'orphan' commas\n", call.=FALSE)

  # check if 'classic' dets file, which has a different column order from the current default
  if(ncol(dets) > 4)
    if(ncol(dets) == 5) { # then probably a 'new' dets file
      if((name[5] %in% cc.names) && min(dets[,5]) >= 0 && max(dets[,5]) <= 4) {} else # extra check for correct values
        stop("unexpected name or values in fifth column (cc, should be between 0 and 4). Please check the manual for guidelines in producing a correct .csv file.\n", call.=FALSE)
    } else
      if(ncol(dets) == 6) { # probably an 'old' file: dR, dSTD, but could also be cc and delta.R (so no column for delta.STD)
        if(name[5] %in% dR.names && name[6] %in% dSTD.names) {
          message("\nHELP!!! 6!!!\n")
          dets <- cbind(dets[,1:4], rep(cc, nrow(dets)), dets[,5:6]) # some shuffling
          message(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD. \nAdding calibration curve column (fifth column, before dR and dSTD) and saving as", csv.file)
          changed <- 1
        } else
	      stop("unexpected names for columns 5/6. If you want to include delta.R, also add a column for delta.STD. Check the manual for guidelines to producing a correct .csv file.\n", call.=FALSE)
      } else
        if(ncol(dets) == 7) { # probably a 'new' file: cc, dR, dSTD
          if(name[5] %in% cc.names && min(dets[,5]) >= 0 && max(dets[,5]) <= 4 &&
            name[6] %in% dR.names && name[7] %in% dSTD.names)
              {} else
                 stop("unexpected column names, order or values in dets file. \nPlease check the manual for correct dets file formats.\n", call.=FALSE)
        } else
          if(ncol(dets) == 8) { # probably an 'old' file: dR, dSTD, ta, tb
            if(name[5] %in% dR.names && name[6] %in% dSTD.names)
            if(name[7] %in% ta.names && name[8] %in% tb.names)
            if(range(dets[,8] - dets[,7]) == c(1,1)) { # check that these set expected student-t values
              dets <- cbind(dets[,1:4], rep(cc, nrow(dets)), dets[,5:6]) # some shuffling
              message(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD. \nAdding calibration curve column (fifth column, before dR and dSTD) and saving as", csv.file)
              changed <- 1
            } else
              stop("unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct .csv file", call.=FALSE)
          } else
            if(ncol(dets) == 9) { # most complex case, many checks needed
              if(name[9] %in% cc.names && # we're almost sure that this is a 'classic' dets file
                min(dets[,9]) >= 0 && max(dets[,9]) <= 4 && # check that this sets calibration curves
                  range(dets[,8] - dets[,7]) == c(1,1) && # check that these set expected student-t values
                    name[5] %in% dR.names && name[6] %in% dSTD.names && # column names as expected?
                      name[7] %in% ta.names && name[8] %in% tb.names) { # column names as expected?
                        dets <- dets[,c(1:4,9,5:8)] # shuffle colums around
                        message(" Assumed order of columns in dets file: lab ID, Age, error, depth, dR, dSTD, t.a, t.b, cc. \nAdapting column order and saving as", csv.file)
                        changed <- 1
                      } else
                        if(name[5] %in% cc.names && # oh, probably a 'new' file from more recent Bacon
                          min(dets[,5]) >= 0 && max(dets[,5]) <= 4 && # check that this sets cal.curves
                            range(dets[,9] - dets[,8]) == c(1,1) && # columns 8-9 set student-t correctly
                              name[8] %in% ta.names && name[9] %in% tb.names && # and are correctly named
                                name[6] %in% dR.names && name[7] %in% dSTD.names) # all lights are green
                                  {} else
                                     stop("unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct .csv file", call.=FALSE)
            } else
              stop("unexpected column names, order or values in dets file. \nPlease check the manual for how to produce a correct dets file.\n", call.=FALSE)

  # more sanity checks
  if(!is.numeric(dets[,2]) || !is.numeric(dets[,3]) || !is.numeric(dets[,4]))
    stop("unexpected values in dets file, I expected numbers. Check the manual.\n", call.=FALSE)
  if(min(dets[,3]) <= 0) {
    message("Warning, zero year errors don't exist in Bacon's world. I will increase them to 1 ", set$age.unit, " yr")
    dets[dets[,3] <= 0,3] <- 1
    changed <- 1
  }
  if( nrow(dets) > 1 && min(diff(dets[,4])) < 0) {
    message("Warning, the depths are not in ascending order, I will correct this")
    dets <- dets[ order(dets[,4]), ]
    changed <- 1
  }

  # if current dets differ from original .csv file, rewrite it
  if(changed > 0)
    write.table(dets, csv.file, sep=paste(sep, "\t", sep=""), dec=dec, row.names=FALSE, col.names=suggested.names[1:ncol(dets)], quote=FALSE)
  dets
}
