#This is a Dockerfile to create an image for shinyproxy
# to build it : sudo docker build -t genepopshiny .
# Add this to the apps entry of application.yml
#  - name: genepop_shiny
#    docker-image: genepopshiny
#    groups: scientists
# Then run : java -jar shinyproxy-0.8.0.jar

FROM openanalytics/r-base

MAINTAINER khalid BELKHIR

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.0.0


# basic R packages
RUN R -e "install.packages(c('devtools', 'git2r'), repos='https://cloud.r-project.org/')"

# basic shiny functionality
RUN R -e "install.packages(c('shiny','shinydashboard','shinyjs', 'rmarkdown','DT', 'yaml'), repos='https://cloud.r-project.org/')"

# install dependencies of the genepop app

RUN R -e "install.packages('ggplot2', repos='https://cloud.r-project.org/')"

#For now the github  you have to git clone "http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop.git" in the current dir
#this contain a copy of genepop-shiny the documentation and sample files
#create a dir to hold Genepop files
#RUN mkdir /root/Genepop
#copy from the host to the image
#COPY Genepop /root/Genepop
#get it from the git repo
RUN cd /root
RUN git clone http://gitlab.mbb.univ-montp2.fr/jlopez/Genepop.git

#install into the docker image
RUN R -e 'devtools::build("/root/Genepop")'
RUN R -e 'devtools::install("/root/Genepop")'

# inside the image move genepop-shiny app and the doc to the /root
RUN mv /root/Genepop/inst/* /root/

COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e shiny::runApp('/root/genepop-shiny/')"]
