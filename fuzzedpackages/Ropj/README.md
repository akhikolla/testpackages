Ropj
====

The goal of this package is to provide the ability to import Origin(R) OPJ
files. The only function, `read.opj(file)`, uses [liborigin] to parse the file
and build a list of its contents. No write support is planned, since it's
absent in [liborigin].

The package is [available on CRAN](https://cran.r-project.org/package=Ropj).

Submodules
----------

If you want to clone the Git repo of this package, don't forget the
`--recursive` flag. Otherwise, use `git submodule update --init --recursive`
after you cloned it.

[liborigin]: https://sourceforge.net/projects/liborigin/
