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
