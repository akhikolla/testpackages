#  fstcore - R bindings to the fstlib library
#
#  Copyright (C) 2017-present, Mark AJ Klik
#
#  This file is part of the fstcore R package.
#
#  The fstcore R package is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Affero General Public License version 3 as
#  published by the Free Software Foundation.
#
#  The fstcore R package is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
#  for more details.
#
#  You should have received a copy of the GNU Affero General Public License along
#  with the fstcore R package. If not, see <http://www.gnu.org/licenses/>.
#
#  You can contact the author at:
#  - fstcore R package source repository : https://github.com/fstpackage/fstcore


#' @useDynLib fstcore, .registration = TRUE
#' @import Rcpp
#' @importFrom utils packageVersion
#' @importFrom utils capture.output
#' @importFrom utils tail
#' @importFrom utils str
#' @importFrom parallel detectCores
NULL


#' R bindings to the fstlib library
#'
#' R package `fstcore` contains R bindings to the C++ `fstlib` library which allows interfacing with _fst_ files.
#' It also contains the `LZ4` and `ZSTD` sources used for compression. `fstcore` exists as a package separate from the
#' `fst` package to facilitate independent updates to the `fstlib`, `LZ4` and `ZSTD` libraries and is used as a
#' backend to `fst`.
#'
#' The fstcore package is based on six C++ libraries or components:
#'
#' * **fstlib**: library containing code to write, read and compute on files stored in the _fst_ format.
#' Written and owned by Mark Klik.
#' * **LZ4**: library containing code to compress data with the LZ4 compressor. Written and owned
#' by Yann Collet.
#' * **ZSTD**: library containing code to compress data with the ZSTD compressor. Written by
#' Yann Collet and owned by Facebook, Inc.
#' * **libdivsufsort-lite**: a  lightweight suffix array construction algorithm. This code is bundled with
#' the ZSTD compression library. Written and owned by Yuta Mori.
#' * **sources from zstdmt**: a multithreading library for Brotli, Lizard, LZ4, LZ5 and Zstandard.  This code is
#' bundled with the ZSTD library. Written and owned by Tino Reichardt.
#' * **source file zstd_opt.h**: Price functions for optimal parser. Written and owned by Yann Collet and Przemyslaw
#' Skibinski. This code is bundled with the ZSTD library.
#'
#' **_The following copyright notice, list of conditions and disclaimer apply to
#' the use of the ZSTD library in the fstcore package:_**
#'
#' BSD License
#'
#' For Zstandard software
#'
#' Copyright (c) 2016-present, Facebook, Inc. All rights reserved.
#'
#' Redistribution and use in source and binary forms, with or without modification,
#' are permitted provided that the following conditions are met:
#'
#' * Redistributions of source code must retain the above copyright notice, this
#' list of conditions and the following disclaimer.
#'
#' * Redistributions in binary form must reproduce the above copyright notice,
#' this list of conditions and the following disclaimer in the documentation
#' and/or other materials provided with the distribution.
#'
#' * Neither the name Facebook nor the names of its contributors may be used to
#' endorse or promote products derived from this software without specific
#' prior written permission.
#'
#' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#' ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#' WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#' DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
#' ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#' (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#' LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
#' ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#' (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#' SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#'
#'
#' **_The following copyright notice, list of conditions and disclaimer apply to
#' the use of the LZ4 library in the fstcore package:_**
#'
#' LZ4 Library
#' Copyright (c) 2011-2016, Yann Collet
#' All rights reserved.
#'
#' Redistribution and use in source and binary forms, with or without modification,
#' are permitted provided that the following conditions are met:
#'
#' * Redistributions of source code must retain the above copyright notice, this
#' list of conditions and the following disclaimer.
#'
#' * Redistributions in binary form must reproduce the above copyright notice, this
#' list of conditions and the following disclaimer in the documentation and/or
#' other materials provided with the distribution.
#'
#' THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#' ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#' WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#' DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
#' ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#' (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#'   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
#' ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#' (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#' SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#'
#'
#' **_The following copyright notice, list of conditions and disclaimer apply to
#' the use of the fstlib library in the fstcore package:_**
#' 
#' fstlib - A C++ library for ultra fast storage and retrieval of datasets
#'
#' Copyright (C) 2017-present, Mark AJ Klik
#'
#' This file is part of fstlib.
#'
#' fstlib is free software: you can redistribute it and/or modify it under the
#' terms of the GNU Affero General Public License version 3 as published by the
#' Free Software Foundation.
#'
#' fstlib is distributed in the hope that it will be useful, but WITHOUT ANY
#' WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#' A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
#' details.
#'
#' You should have received a copy of the GNU Affero General Public License
#' along with fstlib. If not, see <http://www.gnu.org/licenses/>.
#'
#'
#' **_The following copyright notice, list of conditions and disclaimer apply to
#' the use of the libdivsufsort-lite library in the fstcore package:_**
#'
#' Copyright (c) 2003-2008 Yuta Mori All Rights Reserved.
#'
#' Permission is hereby granted, free of charge, to any person
#' obtaining a copy of this software and associated documentation
#' files (the "Software"), to deal in the Software without
#' restriction, including without limitation the rights to use,
#' copy, modify, merge, publish, distribute, sublicense, and/or sell
#' copies of the Software, and to permit persons to whom the
#' Software is furnished to do so, subject to the following
#' conditions:
#'
#' The above copyright notice and this permission notice shall be
#' included in all copies or substantial portions of the Software.
#'
#' THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#' EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#' OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#' NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#' HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#' WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#' FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#' OTHER DEALINGS IN THE SOFTWARE.
#'
#'
#' **_The following copyright notice, list of conditions and disclaimer apply to
#' the use of sources from the zstdmt library included in the fstcore package
#' as part the ZSTD library:_**
#'
#' Copyright (c) 2016 Tino Reichardt
#' All rights reserved.
#'
#' This source code is licensed under both the BSD-style license (found in the
#' LICENSE file in the root directory of this source tree) and the GPLv2 (found
#' in the COPYING file in the root directory of this source tree).
#'
#' You can contact the author at:
#' - zstdmt source repository: https://github.com/mcmilk/zstdmt
#'
#'
#'
#' **_The following copyright notice, list of conditions and disclaimer apply to
#' the use of zstd_opt.h included in the fstcore package as part the ZSTD library:_**
#'
#'
#' Copyright (c) 2016-present, Przemyslaw Skibinski, Yann Collet, Facebook, Inc.
#' All rights reserved.
#'
#' This source code is licensed under the BSD-style license found in the
#' LICENSE file in the root directory of https://github.com/facebook/zstd.
#'
#' This program is free software; you can redistribute it and/or modify it under
#' the terms of the GNU General Public License version 2 as published by the
#' Free Software Foundation. This program is dual-licensed; you may select
#' either version 2 of the GNU General Public License ("GPL") or BSD license
#' ("BSD").


#'
#' @md
#' @docType package
#' @name fstcore-package
#' @aliases fstcore-package
NULL


.onLoad <- function(libname, pkgname) {  # nolint

  # check option for restoring threads after forking
  option_restore <- getOption("fst_restore_after_fork")

  if (!is.null(option_restore)) {
    if (is.logical(option_restore) && !is.na(option_restore)) {
      restore_after_fork(option_restore)
    }
  }

  # set the number of threads here. .onAttach will read threads_fst() and display a
  # useful startup message.
  option_threads <- getOption("fst_threads")

  if (!is.null(option_threads)) {
    if (!is.numeric(option_threads) || is.na(option_threads)) {
      # don't use option if improperly set
      option_threads <- NULL
    }
  }

  if (is.null(option_threads)) {
    logical_cores <- parallel::detectCores(logical = TRUE)

    # if R can't figure out how many logical cores, ask OpenMP to use all
    logical_cores <- ifelse(is.na(logical_cores), 0L, logical_cores)

    # The default number of cores is set to the number of logical cores available on the system.
    # Benchmarks show that hyperthreading increases the read- and write performance of fst.
    setnrofthreads(logical_cores)
  } else {
    # don't need to validate here; threads_fst checks its input
    setnrofthreads(option_threads[1])
  }
}


.onUnload <- function (libpath) {  # nolint
  library.dynam.unload("fstcore", libpath)
}
