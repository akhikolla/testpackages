# MUTEXES

#' @import Rcpp
#' @useDynLib synchronicity

#' @export
setClass('mutex')

#' Lock and Unlock a Mutex
#' 
#' @description The \code{lock} and \code{unlock} functions allow a user to
#' specify exclusive or shared access to a resource.
#' @param m a mutex.
#' @param ... options associated with the mutex being used including
#' \code{block} which forces the mutex to return immediately after trying
#' to acquire a lock
#' @details A call to \code{lock} gives exclusive access to a resource; no other
#' mutex may acquire a lock.  A call to to \code{lock.shared} allows other
#' mutexes to acquire a shared lock on the resource.  When shared lock is
#' called while a exclusive lock has been acquired, the shared lock will
#' block until the exclusive lock is release.  Likewise, if an exclusive lock
#' is called while a shared lock has been acquired, the exclusive lock will
#' block until the shared lock is released.
#' @return The function returns \code{TRUE} if the lock is successfully 
#' called and \code{FALSE} otherwise
#' @examples 
#' m = boost.mutex()
#' lock(m)
#' # Some code that needs to be synchronized...
#' unlock(m)
#' @rdname lock-and-unlock-methods
#' @aliases unlock.shared,boost.mutex-method
#' @export
setGeneric('lock', function(m, ...) standardGeneric('lock'))

#' @rdname lock-and-unlock-methods
#' @export
setGeneric('lock.shared', function(m, ...) standardGeneric('lock.shared'))

#' @rdname lock-and-unlock-methods
#' @export
setGeneric('unlock', function(m, ...) standardGeneric('unlock'))

#' @rdname lock-and-unlock-methods
#' @export
setGeneric('unlock.shared', function(m, ...) standardGeneric('unlock.shared'))

#' @export
setClass('boost.mutex', contains='mutex', 
  representation(isRead='logical', mutexInfoAddr='externalptr'))

#' @title Is it a shared mutex?
#' 
#' @description Tells the user if a mutex is a shared mutex. If it is 
#' not then it must be a write (exclusive) mutex.
#' @docType methods
#' @rdname shared-methods
#' @param m the mutex 
#' @return TRUE if the mutex is shared, FALSE otherwise.
#' @export
setGeneric('shared', function(m) standardGeneric('shared'))

#' @rdname shared-methods
#' @aliases shared,boost.mutex-method
#' @export
setMethod('shared', signature(m='boost.mutex'), function(m) 
  IsRead(m@mutexInfoAddr))

#' @export
setMethod('lock', signature(m='boost.mutex'),
  function(m, ...)
  {
    block = match.call()[['block']]
    if (is.null(block)) block=TRUE
    if (!is.logical(block)) stop('The block argument should be logical')
    if (block) boost_lock(m@mutexInfoAddr)
    else boost_try_lock(m@mutexInfoAddr)
  })

#' @export
setMethod('lock.shared', signature(m='boost.mutex'),
  function(m, ...)
  {
    block = match.call()[['block']]
    if (is.null(block)) block=TRUE
    if (!is.logical(block)) stop('The block argument should be logical')
    if (block) boost_lock_shared(m@mutexInfoAddr)
    else boost_try_lock_shared(m@mutexInfoAddr)
  })

#' @aliases unlock,boost.mutex-method unlock.shared,boost.mutex-method
#' @export
setMethod('unlock', signature(m='boost.mutex'),
  function(m, ...) {
    boost_unlock(m@mutexInfoAddr)
  })

#' @export
setMethod("unlock.shared", signature(m='boost.mutex'),
  function(m, ...) {
    boost_unlock_shared(m@mutexInfoAddr)
  })

#' @export
setGeneric('shared.name', function(m) standardGeneric('shared.name'))

#' @export
setMethod('shared.name', signature(m='boost.mutex'), 
  function(m) GetResourceName(m@mutexInfoAddr) )

#' @export
setGeneric('timeout', function(m) standardGeneric('timeout'))

#' @export
setMethod('timeout', signature(m='boost.mutex'),
  function(m) GetTimeout(m@mutexInfoAddr))

#' @export
setGeneric('is.timed', function(m) standardGeneric('is.timed'))

#' @export
setMethod('is.timed', signature(m='boost.mutex'),
  function(m)
  {
    return(!is.null(timeout(m)))
  })

#' @importFrom uuid UUIDgenerate
#' @export
boost.mutex=function(sharedName=NULL, timeout=NULL, create=TRUE)
{
  isRead = TRUE
  if (!is.null(timeout) && !is.numeric(timeout))
    stop("The timeout parameter must be numeric.")
  if (is.numeric(timeout) && timeout <= 0)
    stop("You must specify a timeout greater than zero.")

  mutexInfoAddr = try({
    if (create) {
      if (is.null(sharedName)) {
        sharedName <- UUIDgenerate()
        # Darwin's shared resources are only 5 characters long.
        if (Sys.info()['sysname'] == "Darwin") {
          sharedName <- substr(sharedName, 1, 5)
        }
      }
      CreateBoostMutexInfo(sharedName, as.double(timeout))
    }
    else
      AttachBoostMutexInfo(sharedName, as.double(timeout))
  })
  
  return(new('boost.mutex', isRead=isRead, mutexInfoAddr=mutexInfoAddr))
}


#' @title An S4 class holding mutex description information.
#' 
#' @description Objects of class description allow users to ``attach'' 
#' to existing mutexes within or across processes.
#' @slot description the list of description information.
#' @export
setClass('descriptor', representation(description='list'))

#' @title Accessor for descriptor objects
#' 
#' @description Retrieve the list of description information from a 
#' descriptor object.
#' @docType methods
#' @rdname description-methods
#' @param x the descriptor object.
#' @return a list of description information.
#' @export
setGeneric('description', function(x) standardGeneric('description'))

#' @rdname description-methods
#' @aliases description,descriptor-method
#' @export
setMethod('description', signature(x='descriptor'),
  function(x) return(x@description))

#' @title An S4 class holding boost.mutex description information.
#' 
#' @description Objects of class description allow users to ``attach'' 
#' to existing mutexes within or across processes.
#' @slot description the list of description information.
#' @export
setClass('boost.mutex.descriptor', contains='descriptor')

#' @title Describe the boost.mutex object
#' 
#' @description The information required to ``attach'' to an existing
#' mutex object.
#' @param x the boost mutex object to describe.
#' @importFrom bigmemory.sri describe
#' @import methods
#' @rdname describe-methods
#' @aliases describe,boost.mutex-method
#' @export
setMethod('describe', signature(x='boost.mutex'),
  function(x)
  {
    return(new('boost.mutex.descriptor',
      description=list(shared.name=shared.name(x),
      timeout=timeout(x))))
  })

#' @title Attach to an existing mutex.
#' 
#' @description Attach to an existing mutex using either a file or 
#' description object
#' @docType methods
#' @rdname attach.mutex-methods
#' @param obj the descriptor object.
#' @param ... other arguments needed by attach.
#' @return A mutex.
#' @export
setGeneric('attach.mutex', function(obj, ...) 
  standardGeneric('attach.mutex'))

#' @rdname attach.mutex-methods
#' @aliases attach.mutex,character-method
#' @export
setMethod('attach.mutex', signature(obj='character'),
  function(obj, ...)
  {
    path = match.call()[['path']]
    if (is.null(path))
    {
      path <- '.'
    }
    path <- path.expand(path)
    if (basename(obj) != obj)
    {

        warning(paste("Two paths were specified in attach.mutex",
          "The one associated with the file will be used.", sep="  "))
      path <- dirname(obj)
      obj <- basename(obj)
    }
    fileWithPath <- file.path(path, obj)
    fi = file.info(fileWithPath)
    print(dir())
    if (is.na(fi$isdir))
      stop( paste("The file", fileWithPath, "could not be found") )
    if (fi$isdir)
      stop( fileWithPath, "is a directory" )
    info <- dget(fileWithPath)
    return(attach.mutex(info, path=path))
  })

#' @rdname attach.mutex-methods
#' @aliases attach.mutex,boost.mutex.descriptor-method
#' @export
setMethod('attach.mutex', signature(obj='boost.mutex.descriptor'),
  function(obj, ...)
  {
    desc = description(obj)
    return(boost.mutex(sharedName = desc$shared.name,
      timeout = desc$timeout, FALSE))
  })


