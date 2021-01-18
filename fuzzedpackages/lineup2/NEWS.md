## lineup2 0.2-5 (2020-10-31)

- Added argument `threshold` to `get_problems()`

- Get rid of rownames from output of `get_problems()`

- Added argument `xlabel` for `hist_self_nonself()` to control x-axis labels.

- Added a vignette

- Added example data

- Added some more examples in the help files.


## lineup2 0.1-7 (2020-10-08)

- Fixed use of `class()` in `is_cluster()`. Instead of
  `"blah" %in% class(object)`, use `inherits(object, "blah")`.


## lineup2 0.1-5 (2019-08-12)

- Revised tests due to change in R's sample() in R 3.6.0.


## lineup2 0.1-4 (2019-03-10)

- Revised documentation to use Markdown throughout.
