Note on installing conda on Ubuntu 24:

mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
source ~/miniconda3/bin/activate
conda init --all
logout
login
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
conda create -y -n mxlib_env
conda activate mxlib_env
bash ./setup/conda_env_setup_Linux.bash

Need to install make, gcc, g++, pkg-config cmake

cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..

On macOS Apple Silicon, if CMake fails its initial compiler test with x86 flags such as `-march=core2` or `-mssse3`,
start from a clean build directory and configure explicitly with the active arm64 conda compiler:

```bash
rm -rf _build
unset CC CXX CFLAGS CXXFLAGS CPPFLAGS LDFLAGS
cmake -S . -B _build \
  -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
  -DCMAKE_C_COMPILER="$(which clang)" \
  -DCMAKE_CXX_COMPILER="$(which clang++)" \
  -DCMAKE_OSX_ARCHITECTURES=arm64
```
