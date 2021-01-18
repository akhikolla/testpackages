FROM r-base:latest

LABEL maintainer "Jean-Francois Rey <jean-francois.rey@inra.fr>"

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev \
    libpq-dev \
    libgeos-dev \
    locales \
    libproj-dev \
    libgdal-dev gdal-bin \
    libudunits2-0 libudunits2-data libudunits2-dev \
    git \
    openssh-client \
    libssh2-1-dev \
    libgit2-dev \
    libglib2.0-dev/unstable \
    libgsl-dev \
    ffmpeg \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

ENV CONTAINER_TIMEZONE Europe/Paris

RUN sudo echo "Europe/Paris" > /etc/timezone
RUN echo "fr_FR.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen fr_FR.UTF8 \
  && /usr/sbin/update-locale LANG=fr_FR.UTF-8

ENV LC_ALL fr_FR.UTF-8
ENV LANG fr_FR.UTF-8

RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" \
  && VERSION=$(cat version.txt) \
  && wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb \
  #&& gdebi -n ss-latest.deb \
  && dpkg -i ss-latest.deb \
  && rm -f version.txt ss-latest.deb \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('shiny','htmltools', 'DT', 'shinyjs', 'Rcpp', 'rgdal', 'sp', 'devtools', 'Matrix', 'mvtnorm', 'rgeos', 'maptools', 'fields', 'splancs', 'sf', 'RSQLite', 'foreach', 'parallel', 'doParallel', 'gridExtra', 'png', 'grid', 'future', 'promises'), repos='https://cran.biotools.fr/')"

RUN echo 'run_as shiny;\n\
  preserve_logs true;\n\
  server { \n\
  listen 3838;\n\
  location / { \n\
    app_init_timeout 120; \n\
    app_idle_timeout 30; \n\
    site_dir /srv/shiny-server/landsepi; \n\
    log_dir /var/log/shiny-server;\n\
    directory_index off; \n\
  }\n\
}' > /etc/shiny-server/shiny-server.conf

RUN echo '#!/bin/sh\n\
  touch /var/log/shiny-server.log\n\
  chown shiny.root /var/log/shiny-server.log\n\
  touch /var/run/shiny-server.pid\n\
  chown shiny.root /var/run/shiny-server.pid\n\
  chown -R shiny.shiny /var/lib/shiny-server\n\
  su shiny -c"shiny-server --pidfile=/var/run/shiny-server.pid >> /var/log/shiny-server.log 2>&1"\n' > /usr/bin/shiny-server.sh

COPY . /landsepi/
RUN R CMD build --no-build-vignettes /landsepi \
  && R CMD INSTALL --no-docs landsepi_*.tar.gz \
  && rm -r /landsepi

RUN apt-get update && apt-get remove --purge -y \
    gdebi-core \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev \
    libpq-dev \
    libgeos-dev \
    libproj-dev \
    libgdal-dev \
    libudunits2-dev \
    git \
    openssh-client \
    libssh2-1-dev \
    libgit2-dev \
  && apt-get autoremove -y \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

#USER shiny
RUN passwd shiny -d
RUN mkdir -p /var/log/shiny-server \
  && chown shiny.shiny /var/log/shiny-server

RUN mkdir -p /srv/shiny-server/landsepi
COPY inst/shiny-landsepi/ui.R /srv/shiny-server/landsepi
COPY inst/shiny-landsepi/server.R /srv/shiny-server/landsepi
COPY inst/shiny-landsepi/global.R /srv/shiny-server/landsepi
COPY inst/shiny-landsepi/www/* /srv/shiny-server/landsepi/www/
RUN chown -R shiny.root /srv/shiny-server
RUN chown -R shiny.root /usr/bin/shiny-server.sh \
  && chmod u+x /usr/bin/shiny-server.sh


WORKDIR /home/shiny

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]


