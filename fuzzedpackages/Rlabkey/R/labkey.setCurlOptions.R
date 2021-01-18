##
#  Copyright (c) 2014-2018 LabKey Corporation
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
##
PACKAGE_ENV = new.env()

# Helper function to get the base curl options used for all http or https requests
#
labkey.setCurlOptions <- function(...)
{
    # test for legacy config names
    params <- c(...)
    if (is.element('ssl.verifyhost', names(params)))
    {
        stop(paste("The legacy config : ssl.verifyhost is no longer supported please update to use : ssl_verifyhost"))
    }

    if (is.element('ssl.verifypeer', names(params)))
    {
        stop(paste("The legacy config : ssl.verifypeer is no longer supported please update to use : ssl_verifypeer"))
    }

    # default curl options
    options <- config(ssl_verifyhost=2, ssl_verifypeer=TRUE, followlocation=TRUE, sslversion=1L, useragent = "Rlabkey")

    # check if a certificate bundle has been specified from the environment variable
    vars <- Sys.getenv("RLABKEY_CAINFO_FILE")
    if (nchar(vars[1]) > 0)
    {
        options <- c(options, config(cainfo = vars[1]))
    }

    # merge in any overrides
    options <- c(options, config(...))
    assign("RLABKEY_CURL_OPTIONS", options, envir=PACKAGE_ENV)

    return(get("RLABKEY_CURL_OPTIONS", envir=PACKAGE_ENV))
}

labkey.acceptSelfSignedCerts <- function()
{
    return(labkey.setCurlOptions(ssl_verifyhost=0, ssl_verifypeer=FALSE))
}
