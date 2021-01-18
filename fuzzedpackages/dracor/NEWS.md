# dracor 0.2.4

* fix compilation on gcc11 due to missing #include in draco library
  (reported by Uwe Ligges when 0.2.3 was submitted)

# dracor 0.2.3

* fix compilation on gcc11 due to missing #include in draco library (#3)
  (reported by BDR)
* fix for README URL
* dev: travis fixes

# dracor 0.2.2

* remove excess Draco library source code files (#1)
* fix build error on windows-oldrel due to gcc 4.9.3 compile error (#2)

# dracor 0.2.1

* `draco_decode()` can now decode a URL or file on disk
* Add `Enhances: rgl` so that we can use rgl in examples
* Better error messages

# dracor 0.2

* First released version
* Passes local `R CMD check` without notes
* default return object is `rgl::mesh3d` (but package does not depend on rgl)
* now includes basic README

# dracor 0.1

* Added a `NEWS.md` file to track changes to the package.
