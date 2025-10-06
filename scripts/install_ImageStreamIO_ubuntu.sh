#!/bin/bash
set -o pipefail
# since we care about almost every exit status, check them all. use `|| true` to bypass.
set -e
#################################################################################
# A script to install ImageStreamIO on Ubuntu 24.04 (and probably other versions)
#
#
# Creates a directory called ~Source TODO: make this configurable
# Switches ISIO repo to dev TODO: make this configurable
#################################################################################

sudo apt-get install -y \
     git \
     make \
     cmake \
     pkg-config \
     gcc \
     g++ \
     libcfitsio-dev


## Make work directory
mkdir -p $HOME/Source
cd $HOME/Source

## Install ImageStreamIO
if [[ ! -d ImageStreamIO ]]; then
    git clone https://github.com/milk-org/ImageStreamIO.git

    cd ImageStreamIO
    git checkout dev
    mkdir -p _build
    cd _build
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
    make
    sudo -E make install
fi

