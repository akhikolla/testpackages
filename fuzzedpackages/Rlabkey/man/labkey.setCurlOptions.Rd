\name{labkey.setCurlOptions}
\alias{labkey.setCurlOptions}
\title{Modify the current set of Curl options that are being used in the existing session}
\description{
Rlabkey uses the package httr to connect to the LabKey Server.}
\arguments{
  \item{options args}{a variable list of arguments to set the RCurl options}
  \item{ssl_verifyhost}{check the existence of a common name and also verify that it matches the hostname provided}
  \item{ssl_verifypeer}{specifies whether curl will verify the peer's certificate}
  \item{followlocation}{specify is curl should follow any location header that is sent in the HTTP request}
  \item{sslversion}{the SSL version to use}
}

\details{
This topic explains how to configure Rlabkey to work with a LabKey Server running SSL.\cr

Rlabkey uses the package httr to connect to the LabKey Server. On Windows, the httr package is not configured for SSL
by default. In order to connect to a HTTPS enabled LabKey Server, you will need to perform the following steps:\cr

1. Create or download a "ca-bundle" file.\cr

We recommend using ca-bundle file that is published by Mozilla. See http://curl.haxx.se/docs/caextract.html. You have two options:\cr

Download the ca-bundle.crt file from the link named "HTTPS from github:" on http://curl.haxx.se/docs/caextract.html\cr
Create your own ca-bundle.crt file using the instructions provided on http://curl.haxx.se/docs/caextract.html\cr

2. Copy the ca-bundle.crt file to a location on your hard-drive.\cr
If you will be the only person using the Rlabkey package on your computer, we recommend that you\cr

create a directory named `labkey` in your home directory\cr
copy the ca-bundle.crt into the `labkey` directory\cr

If you are installing this file on a server where multiple users will use may use the Rlabkey package, we recommend that you
create a directory named `c:labkey`\cr

copy the ca-bundle.crt into the `c:labkey` directory\cr

3. Create a new Environment variable named `RLABKEY_CAINFO_FILE`\cr

On Windows 7, Windows Server 2008 and earlier\cr

Select Computer from the Start menu.\cr
Choose System Properties from the context menu.\cr
Click Advanced system settings > Advanced tab.\cr
Click on Environment Variables.\cr
Under System Variables click on the new button.\cr
For Variable Name: enter RLABKEY_CAINFO_FILE\cr
For Variable Value: enter the path of the ca-bundle.crt you created above.\cr
Hit the Ok buttons to close all the windows.\cr
On Windows 8, Windows 2012 and above\cr
Drag the Mouse pointer to the Right bottom corner of the screen.\cr
Click on the Search icon and type: Control Panel.\cr
Click on -> Control Panel -> System and Security.\cr
Click on System -> Advanced system settings > Advanced tab.\cr
In the System Properties Window, click on Environment Variables.\cr
Under System Variables click on the new button.\cr
For Variable Name: enter RLABKEY_CAINFO_FILE\cr
For Variable Value: enter the path of the ca-bundle.crt you created above.\cr
Hit the Ok buttons to close all the windows.\cr
Now you can start R and begin working.
}
