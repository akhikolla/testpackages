## Internal errors ----

stop_bad_class = function(x, obj_name) {
    stop("`x` must be an `", obj_name,"` object. Class ", class(x)[1], " is not supported.")
}
