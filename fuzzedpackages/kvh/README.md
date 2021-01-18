# kvh
This package can read/write KVH files in R.
The format [KVH](http://serguei.sokol.free.fr/kvh-format/) is a lightweight format that can be read/written both by humans and machines.
It can be useful in situations where XML or alike formats seem to be an overkill.
We provide an ability to parse KVH files in R pretty fast due to 'Rcpp' use.
Key/Values that are returned by kvh_read() are always character strings. User has to convert them furthermore to somewhat usefull for him.

Example:
```
     # prepare object to write to kvh file
     obj=list(x=structure(1:3, names=letters[1:3]), R=R.version)
     # write it
     obj2kvh(obj, "test", "test.kvh") # will create test.kvh file
     # read it back
     l=kvh_read("test.kvh")
     # check a field
     l$test$x # NB. it has a character values put in a list not a numeric vector as it was in obj.
     attr(l$test$x, "ln") # line number where the entry test/x started in test.kvh
```
