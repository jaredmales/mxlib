#!/bin/bash
set -o pipefail
# since we care about almost every exit status, check them all. use `|| true` to bypass.
set -e
#################################################################################
# A script to provision an Ubuntu 24.04 machine with mxlib and its dependencies
# Note: this probably works with 22.04 too
#
#
# Creates a directory called ~Source TODO: make this configurable
# Switches some repos to dev TODO: make this configurable
#################################################################################

sudo apt-get install -y \
    git \
    make \
    cmake \
    pkg-config \
    gcc \
    g++ \
    libgsl-dev \
    libboost-all-dev \
    libcfitsio-dev \
    libopenblas-dev \
    liblapacke-dev \
    libfftw3-dev \
    libeigen3-dev

## Make work directory
mkdir -p $HOME/Source


## mxlib
cd $HOME/Source
if [[ ! -d mxlib ]]; then
    git clone https://github.com/jaredmales/mxlib.git
fi
cd mxlib
git checkout dev
mkdir -p _build
cd _build
cmake ..
make
sudo make install







