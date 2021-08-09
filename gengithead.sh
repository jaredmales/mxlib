#!/bin/bash

##############################################################################
# gengithead.sh: a script producing a .h file with git repo information
#
# Copyright (C) 2014-2021 Jared R. Males
#
# This script is licensed under the terms of the MIT license.
# https://opensource.org/licenses/MIT
#
# Author: Jared Males (jaredmales@pm.me)
# 
# Contributors: Joseph Long
#
# Description: This script produces a C/C++ header (.h) with information
#              about the state of a git repository. This allows recording
#              of the state of a git repo at compilation, including whether
#              or not the repo was modified.
#            
#              This also prints a message at compile time to warn if the repo
#              was modified.
#
#              The goal is allow for reproducibility.  For instance, I write
#              the git hash to FITS image headers processed by my code.
#
# Usage: gengithead.sh directory output_path prefix
#
#        directory = directory of git repo, optional default is './'
#        output_path = optional ouput path, default './git_version.h'
#        prefix = optional prefix for #define variable names, default is GIT
#
#        note: the arguments are positional, so to change output_path you must 
#              define directory first.
##############################################################################

set -euo pipefail

#center text in a the pragma message argument to be pretty
centerName(){
  textsize=${#1}
  width=42
  span=$(((width + textsize) / 2 ))
  espan=$((span - textsize))
  if [ $((espan % 2)) -eq 1 ]; then espan=$((espan - 1)); fi
  printf "    #pragma message (\"*%${span}s%${espan}s\")\n" "$1" "*"
}

#defaults for args
GITPATH=${1:-'./'}
HEADPATH=${2:-'./git_version.h'}
PREFIX=${3:-'GIT'}

GIT_HEADER="$HEADPATH"

GIT_BRANCH=$(git --git-dir=$GITPATH/.git --work-tree=$GITPATH rev-parse --abbrev-ref HEAD)
GIT_VERSION=$(git --git-dir=$GITPATH/.git --work-tree=$GITPATH log -1 --format=%H)

set +e
git --git-dir=$GITPATH/.git --work-tree=$GITPATH diff-index --quiet HEAD --
GIT_MODIFIED=$?
set -e

REPO_NAME=$(basename $(git --exec-path=$GITPATH rev-parse --show-toplevel))

echo "#ifndef $PREFIX""_VERSION_H" > $GIT_HEADER
echo "#define $PREFIX""_VERSION_H" >> $GIT_HEADER
echo "" >> $GIT_HEADER
echo "#define $PREFIX""_BRANCH \"$GIT_BRANCH\"" >> $GIT_HEADER
echo "#define $PREFIX""_CURRENT_SHA1 \"$GIT_VERSION\"" >> $GIT_HEADER
echo "#define $PREFIX""_REPO_MODIFIED  $GIT_MODIFIED"  >> $GIT_HEADER
echo "" >> $GIT_HEADER
echo "" >> $GIT_HEADER
if [ $GIT_MODIFIED = 1 ]; then
echo "#if $PREFIX""_REPO_MODIFIED == 1" >> $GIT_HEADER
echo "  #ifndef GITHEAD_NOWARNING" >> $GIT_HEADER
echo "    #pragma message (\"******************************************\")" >> $GIT_HEADER
echo "    #pragma message (\"*                                        *\")" >> $GIT_HEADER
centerName "WARNING: repository modified" >> $GIT_HEADER
centerName "changes not committed for"  >> $GIT_HEADER
centerName $REPO_NAME >> $GIT_HEADER
echo "    #pragma message (\"*                                        *\")" >> $GIT_HEADER
echo "    #pragma message (\"******************************************\")" >> $GIT_HEADER
echo "  #endif" >> $GIT_HEADER
echo "#endif" >> $GIT_HEADER
fi
echo "" >> $GIT_HEADER
echo "" >> $GIT_HEADER
echo "#endif" >> $GIT_HEADER
