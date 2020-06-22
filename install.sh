#!/bin/bash

#ensure pip is installed and runnable from command line
pip --version &> /dev/null || {
    echo "Error: pip python package manager not found.  Install pip using your package manager to continue. Exiting now." 
    exit
}

#ensure G++ compiler is installed
g++ --version &> /dev/null || {
    echo "Error: G++ compiler not found.  Install the GNU C++ compiler using your package manager to continue. Exiting now." 
    exit
}

ROOT_DIR=`pwd`
PATH_EXPORT_CMD="export PATH=\$PATH:$ROOT_DIR/bin/"

#if an argument is given, use it as the login profile 
[ $# -ne 0 ] && PROFILE_FILE=$1 || PROFILE_FILE="$HOME/.profile" 

#If not already added, append SPHERLSanal to path in profile file
grep -q $PATH_EXPORT_CMD $PROFILE_FILE || sed -i.bak "$ a $PATH_EXPORT_CMD" $PROFILE_FILE 

source $PROFILE_FILE

#install python packages using pip
pip install h5py
pip install pyevtk
pip install eos-py
pip install numpy

#build modules from source code
cd ./src
make
make clean

#Make working directories
mkdir eos
mkdir data