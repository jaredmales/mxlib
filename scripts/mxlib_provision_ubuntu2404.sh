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

# from https://github.com/milk-org/milk/blob/dev/Dockerfile
sudo apt-get install -y \
    git \
    make \
    dpkg-dev \
    libc6-dev \
    cmake \
    pkg-config \
    python3-dev \
    libcfitsio-dev \
    pybind11-dev \
    python3-pybind11 \
    libgsl-dev \
    libfftw3-dev \
    libncurses-dev \
    libbison-dev \
    libfl-dev \
    libreadline-dev \
    gfortran libopenblas-dev liblapacke-dev \
    pkg-config \
    gcc \
    g++


## Make work directory
mkdir -p ~/Source
cd ~/Source

## Install milk
if ! command -v milk; then
    if [[ ! -d milk ]]; then
        git clone --recursive https://github.com/milk-org/milk.git
    fi
    cd milk
    git checkout dev
    mkdir -p _build
    cd _build
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
    sudo make install
    if [[ ! -h /usr/local/milk ]]; then
        sudo ln -sfv /usr/local/milk-* /usr/local/milk
    fi
fi

## Setup milk for linking
if [[ ! -z /etc/ld.so.conf.d/milk.conf ]]; then
    echo /usr/local/milk/lib/ | sudo tee -a  /etc/ld.so.conf.d/milk.conf
    sudo ldconfig
fi

## mxlib
cd ~/Source
sudo apt-get install -y libgsl-dev libboost-all-dev libcfitsio-dev libopenblas-dev libfftw3-dev libeigen3-dev
if [[ ! -d mxlib ]]; then
    git clone https://github.com/jaredmales/mxlib.git
fi
cd mxlib
git checkout dev
if [[ ! -z local/Common.mk ]]; then
    echo NEED_CUDA=no > local/Common.mk
    echo PREFIX=/usr/local >> local/Common.mk
fi
make
sudo make install

