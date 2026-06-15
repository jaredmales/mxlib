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
    pkg-config \
    gcc \
    g++ \
    libgsl-dev \
    libboost-system-dev \
    libcfitsio-dev \
    libopenblas-dev \
    liblapacke-dev \
    libfftw3-dev \
    libeigen3-dev

#make sure old cmake is gone and upgrade to latest
sudo apt purge cmake
sudo snap install cmake --classic
hash -r

#get an up-to-date g++ with format
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install gcc-13 g++-13
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 130 --slave /usr/bin/g++ g++ /usr/bin/g++-13   

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







