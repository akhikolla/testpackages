
.onAttach = function(libname, pkgname)
{
    temp = packageDescription("lazygreedy")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version,
                " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "copyright (c) 2020, Bokgyeong Kang, John Hughes, Quirijn W. Bouts, Alex P. ten Brink, and Kevin Buchin\n",
                sep = "")
    msg = paste(msg, 'For citation information, type citation("lazygreedy").\n', sep = "")
    msg = paste(msg, 'Type help(package = lazygreedy) to get started.\n', sep = "")
    packageStartupMessage(msg)
}
