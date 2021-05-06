#!/usr/bin/env bash

function install_common_packages() {
    conda install -y pkg-config cmake cfitsio cppzmq eigen boost mkl mkl-include gsl fftw
    conda env config vars set MKLROOT="$CONDA_PREFIX"
}
