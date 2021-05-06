#!/usr/bin/env bash
set -euo pipefail
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/conda_env_setup.bash
install_common_packages
conda install clang clangxx clang_osx-64 clangxx_osx-64
conda env config vars set CXX=clang++
