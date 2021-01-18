## ---- eval = FALSE------------------------------------------------------------
#  useUnicode(FALSE)

## ---- include=FALSE-----------------------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
library("set6")
set.seed(42)
useUnicode(TRUE)

## -----------------------------------------------------------------------------
empty = Set$new()
empty
summary(empty)
empty$properties
empty$traits

## -----------------------------------------------------------------------------
Set$new(1,2,3)
Set$new(letters[1:5])
Set$new(1, 2i, "a", Set$new(1))

## -----------------------------------------------------------------------------
# A Set cannot have duplicated elements, and ordering does not matter
Set$new(1,2,2,3) == Set$new(3,2,1)

# A Tuple can have duplicated elements, and ordering does matter
Tuple$new(1,2,2,3) != Tuple$new(1,2,3)
Tuple$new(1,3) != Tuple$new(3,1)

# An interval can be an interval of integers or numerics, but must be continuous
Interval$new(1, 10, class = "integer")
Interval$new(1, 10) # numeric is default
# `type` is used to specify the interval upper and lower closure
Interval$new(1, 10, type = "()") 

# SpecialSets are useful for common 'special' mathematical sets
# Use listSpecialSets() to see which are available.
Reals$new()
PosIntegers$new()

# ConditionalSets are user-built sets from logical statements.
# For example, the set of even numbers.
ConditionalSet$new(function(x) x %% 2 == 0)

# Finally FuzzySets and FuzzyTuples expand Sets and Tuples to allow for partial 
# membership of elements. These have two constructors, either by specifying the elements
# and membership alternatively, or by passing these separately to the given arguments.
FuzzySet$new(1, 0.1, 2, 0.2, "a", 0.3) ==
  FuzzySet$new(elements = c(1,2,"a"), membership = c(0.1,0.2,0.3))

## -----------------------------------------------------------------------------
s = Set$new(1,2,3)
s$contains(1)
s$contains(2, 4)
c(2, 4) %inset% s

s$isSubset(Set$new(1,2,3), proper = FALSE)
s$isSubset(Set$new(1,2,3), proper = TRUE)
c(Set$new(1), Set$new(4, 5)) < s

# Sets are FuzzySets with membership = 1
s$equals(FuzzySet$new(elements = 1:3, membership = 1))
s$equals(FuzzySet$new(elements = 1:3, membership = 0.1))
s == Set$new(1, 2, 3)
s != c(Set$new(1,2,3), Set$new(1, 2))

1:10 %inset% ConditionalSet$new(function(x) x %% 2 == 0)

# The `bound` argument in `isSubset` is used for specifying
#   how open interval containedness should be checked
i = Interval$new(1, 10, type = "(]")
i$contains(Set$new(1), bound = FALSE)
i$contains(Set$new(10), bound = FALSE)
i$contains(Set$new(1), bound = TRUE)
i$contains(Set$new(10), bound = TRUE)

## -----------------------------------------------------------------------------
Set$new(1) + Set$new(2) + Set$new(3)
Interval$new(1, 10) + Set$new(1) # no effect
setunion(Set$new(1,2), Interval$new(3, 10), Set$new(16))

PosReals$new() + NegReals$new()

## -----------------------------------------------------------------------------
Set$new(elements = 1:10) - Set$new(elements = 4:10)
Set$new(1,2,3,4) - Set$new(2)

Reals$new() - PosReals$new()

Interval$new(5, 10) - Interval$new(3, 12)
Interval$new(5, 10) - Interval$new(7, 12)
Interval$new(5, 10) - Interval$new(11, 12) # no effect

## -----------------------------------------------------------------------------
Set$new(1, 2) * Set$new(3, 4)
Set$new(1, 2) * Set$new(3, 4) * Set$new(5, 6) # n-ary

# nest = FALSE default - we will return to the `simplify` argument below
setproduct(Set$new(1, 2), Set$new(3, 4), Set$new(5, 6), nest = TRUE, simplify = TRUE)
setproduct(Set$new(1, 2), Set$new(3, 4), Set$new(5, 6), nest = FALSE, simplify = TRUE)

# n-ary cartesian product on the same set
setpower(Set$new(1,2), 3, simplify = TRUE)

## -----------------------------------------------------------------------------
Set$new(1,2,3) & Set$new(3,5,6)
Set$new(5,6) & Reals$new()
Set$new(1) & Set$new(2)

## -----------------------------------------------------------------------------
# default: simplify = TRUE
setunion(Set$new(1,2,3), Set$new(4,5))
setunion(Set$new(1,2,3), Set$new(4,5), simplify = FALSE)

# default: simplify = FALSE
setproduct(Set$new(1,2), Set$new(4,5))
setproduct(Set$new(1,2), Set$new(4,5), simplify = TRUE)

# default: simplify = FALSE
powerset(Set$new(1,2,3))
powerset(Set$new(1,2,3), simplify = TRUE)

## -----------------------------------------------------------------------------
u = setunion(Set$new(1,2,3), Set$new(4,5), simplify = FALSE)
c(2,5,8) %inset% u

p = Set$new(1,2) * Set$new(3,4)
p$contains(Tuple$new(2, 4))

