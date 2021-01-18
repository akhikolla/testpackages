# mnis 0.2.7

## Bug Fixes

Fixed bug in `mnis_tidy()` that didn't correct data schema names.

Fixed test errors for `mnis_joined_between()`.

Removed `basic_details` parameter from `mnis_extra()`, as it didn't provide any information that wasn't included by other parameters.

Changed `mnis_extra()` code to eliminate duplicate columns.

# mnis 0.2.6

Now accepts dates as character values in "YYYY-MM-DD" format, and objects of class Date, POSIXt, POSIXct, POSIXlt or anything else than can be coerced to a date with `as.Date()`.

Fixed bugs on `mnis_joined_between` to return full, proper tibble.

Dropped parameter 'joined_since' from `mnis_all_members`, in favour of `mnis_joined_between`.

Tidied up code in `mnis_reference` functions.

Fixed bug in `mnis_party_state` which produced an error if not using the current system date.

Added optional 'tidy_style' parameter, allowing users to decide which style `"snake_case"`, `"camelCase"` and `"period.case"` they want variable names to be in, if 'tidy'==TRUE.

# mnis 0.2.5

`mnis_mps_on_date` and `mnis_peers_on_date` functions introduced, which return all members of the House of Commons or the House of Lords eligible to sit on a given date.

# mnis 0.2.4

## bug fixes

Bug fixes for constituency names in `mnis_eligible`

# mnis 0.2.3

## tibbling

`mnis` now uses tibbles instead of data frames as the data class returned from API calls.

Bug fixes for constituency names in `mnis_eligible`

# mnis 0.2.2

## Removing BOM, changing behaviour for empty mnis_additional function calls

All functions now remove byte order marks from the API response.

All empty `mnis_additional` functions now default to returning `mnis_all_members` instead of an error.

The `mnis_all_members` without any parameters now returns an error message that tells the user to include data for the ID parameter, rather than a generic API error.


# mnis 0.2.1

## Bug fixes, further rollout of tidy function

Expansion of `tidy` function to all other package functions and improved speed and quality.

Fixed errors in `mnis_extra` that prevented it from functioning properly.

Improved handling of deprecated functions.

New `mnis_all_members` function that returns all members from both houses, with various options accepted as function parameters.


# mnis 0.2.0

## Introducing New Naming Conventions, mnis_eligible and tidy

The old function names were awkward. The new ones are all lower case and use underscores to separate individual words. All old function names have been deprecated, but can still be called.

`mnis_eligible` returns a data frame with information on all members eligible or previously eligible to sit in the House of Lords, the House of Commons, or both.

The Members' Names Information Service returns variables names with extra text, including periods, @ signs and superfluous text stuck on the end of names. The `tidy` parameter gives an option to remove this extra text.


# mnis 0.1.1

## Introducing mnis_JoinedBetween

`mnis_JoinedBetween` returns a data frame with information on all members who joined the House of Lords, the House of Commons, or both, between two given dates.

As the Members' Names Information Service has dozens of different search parameters and hundreds of possible combinations of searches, it is not reasonable to build functions for every possible query to the API. However, I will be rolling out functions for potentially common queries as and when I can. If you have any particular queries you would like functions for please let me know and I'll create it.


# mnis 0.1.0

## Introducing the mnis package

This is the first release of the `mnis` package, which provides functions to download data from the Members' Names Information Service for the UK Houses of Parliament.

See the package documentation for details on each function and the type of data available.
