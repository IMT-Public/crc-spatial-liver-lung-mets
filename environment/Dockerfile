# Code Ocean generates this file from the Environment editor. It is checked in
# so that the capsule can also be rebuilt outside Code Ocean.
#
# Base image matches the R version reported in the manuscript Methods (4.5.1).
FROM rocker/r-ver:4.5.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
      libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
      libtiff5-dev libjpeg-dev libhdf5-dev libgit2-dev zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

COPY install.R /tmp/install.R
RUN Rscript /tmp/install.R

