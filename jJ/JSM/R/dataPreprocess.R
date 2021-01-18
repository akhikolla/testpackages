# Data Preprocess

dataPreprocess <- function(long, surv, id.col, long.time.col, surv.time.col, surv.event.col, 
                           surv.event.indicator = list(censored = 0, event = 1), suffix = ".join") {
  call <- match.call()
  
  if (!(inherits(long, "data.frame") && inherits(surv, "data.frame"))) {
    stop("\n 'long' and 'surv' must be of class 'data.frame'.")
  }
  
  if (!inherits(id.col, "character")) {
    stop("\n 'id.col' must be of type 'character'.")
  }
  
  if (!inherits(long.time.col, "character")) {
    stop("\n 'long.time.col' must be of type 'character'.")
  }
  
  if (!inherits(surv.time.col, "character")) {
    stop("\n 'surv.time.col' must be of type 'character'.")
  }
  
  if (!inherits(surv.event.col, "character")) {
    stop("\n 'surv.event.col' must be of type 'character'.")
  }
  
  if (!(id.col %in% names(long) && id.col %in% names(surv)) || length(id.col) != 1) {
    stop("\n 'id.col' must correspond to the subject identification column in both 'long' and 'surv'.")
  }
  
  if (!(long.time.col %in% names(long)) || length(long.time.col) != 1) {
    stop("\n 'long.time.col' must correspond to the column specifying the time of measurment in 'long'.")
  }
  
  if (!(surv.time.col %in% names(surv)) || length(surv.time.col) != 1) {
    stop("\n 'surv.time.col' must correspond to the column specifying the time of event in 'surv'.")
  }
  
  if (!(surv.event.col %in% names(surv)) || length(surv.event.col) != 1) {
    stop("\n 'surv.event.col' must correspond to the column specifying wheter there is an event in 'surv'.")
  }
  
  tempEvent = unique(surv[[surv.event.col]])
  if (length(tempEvent) > 2) {
    stop("\n 'surv.event.col' column in 'surv' should have no more than two distinct values (censored v.s. event).")
  } else if (!all(tempEvent %in% c(surv.event.indicator$censored, surv.event.indicator$event))) {
    stop("\n values in the 'surv.event.col' column of 'surv' does not match those specified in 'surv.event.indicator.")
  }
  
  tempData <- merge(long, surv, by = 'ID', suffixes = c(".long", ".surv"))
  if (!long.time.col %in% names(tempData)) {
    long.time.col <- paste(long.time.col, '.long', sep = '')
  }
  
  if (!surv.time.col %in% names(tempData)) {
    surv.time.col <- paste(surv.time.col, '.surv', sep = '')
  }
  
  if (!surv.event.col %in% names(tempData)) {
    surv.event.col <- paste(surv.event.col, '.surv', sep = '')
  }
  
  tempData[paste('start', suffix, sep = '')] <- tempData[long.time.col]
  
  ID <- as.vector(tempData[[id.col]])
  tempID <- which(!duplicated(ID))
  tempID <- c(tempID, length(ID) + 1)
  ni <- diff(tempID) 
  end <- c()
  event <- c()
  for (i in 1:length(ni)) {
    if (ni[i] == 1) {
      tempEnd <- tempData[[surv.time.col]][tempID[i]]
      tempEvent <- tempData[[surv.event.col]][tempID[i]]
    } else {
      tempEnd <- c(tempData[[long.time.col]][(tempID[i] + 1):(tempID[i + 1] - 1)], tempData[[surv.time.col]][tempID[i + 1] - 1])
      tempEvent <- c(rep(surv.event.indicator$censored, ni[i] - 1), tempData[[surv.event.col]][tempID[i + 1] - 1])
    }
    end <- c(end, tempEnd)
    event <- c(event, tempEvent)
  }
  
  tempData[paste('stop', suffix, sep = '')] <- end
  tempData[paste('event', suffix, sep = '')] <- event
  
  return(tempData)
}
