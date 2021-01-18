stop.redef <-
function (locstring = "", ...)
{
    print(locstring, quote = FALSE)
    stop(locstring, ...)
}
