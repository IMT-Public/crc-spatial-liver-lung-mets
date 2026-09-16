# Code Ocean generates this file from the Environment editor. It is checked in
# so that the capsule can also be rebuilt outside Code Ocean.
#
# Code Ocean only accepts its own base images. This one ships R 4.5.1, the
# version reported in the manuscript Methods.
FROM registry.codeocean.com/codeocean/r-studio:r-4.5.1-1.2204.0-rstudio-server-2025.05.1-513-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
      libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
      libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
      libtiff5-dev libjpeg-dev libhdf5-dev libgit2-dev zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

COPY install.R /tmp/install.R
RUN Rscript /tmp/install.R
