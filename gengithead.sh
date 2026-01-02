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
#        output_path = optional output path, default './git_version.h'
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
  espan=$((span - textsize+2))
  if [ $((espan % 2)) -eq 1 ]; then espan=$((espan - 1)); fi
  printf "                     \"*\033[33m%${span}s\033[39m%${espan}s\"%s\n" "$1" "*\n" "\\"
}

#defaults for args
GITPATH=${1:-'./'}
REPO_NAME=$(basename $(git --exec-path=$GITPATH rev-parse --show-toplevel))
REPO_PATH=$(realpath $GITPATH)

#default PREFIX is the repo name
PREFIX=${3:-$REPO_NAME}

#default output path is PREFIX_git_version.h
HEADPATH=${2:-"./"$PREFIX"_git_version.h"}
GIT_HEADER="$HEADPATH"

# Get git status
GIT_URL=$(git --git-dir=$GITPATH/.git --work-tree=$GITPATH remote get-url origin)
GIT_BRANCH=$(git --git-dir=$GITPATH/.git --work-tree=$GITPATH rev-parse --abbrev-ref HEAD)
GIT_VERSION=$(git --git-dir=$GITPATH/.git --work-tree=$GITPATH log -1 --format=%H)

# Check if repo is modified
set +e
git --git-dir=$GITPATH/.git --work-tree=$GITPATH diff-index --quiet HEAD --
GIT_MODIFIED=$?
set -e

# Check if untracked files exist
set +e
GIT_UNTRACKED_STR=$(git status | grep Untracked)
set -e

GIT_UNTRACKED=0
if [ ${#GIT_UNTRACKED_STR} -gt 0 ]; then
    GIT_UNTRACKED=1
fi

echo "#ifndef $PREFIX""_GIT_VERSION_H" > $GIT_HEADER
echo "#define $PREFIX""_GIT_VERSION_H" >> $GIT_HEADER
echo "" >> $GIT_HEADER
echo "#define $PREFIX""_REPO \"$REPO_NAME\"" >> $GIT_HEADER
echo "#define $PREFIX""_URL \"$GIT_URL\"" >> $GIT_HEADER
echo "#define $PREFIX""_BRANCH \"$GIT_BRANCH\"" >> $GIT_HEADER
echo "#define $PREFIX""_SRCPATH \"$REPO_PATH\"" >> $GIT_HEADER
echo "#define $PREFIX""_CURRENT_SHA1 \"$GIT_VERSION\"" >> $GIT_HEADER
echo "#define $PREFIX""_REPO_MODIFIED  $GIT_MODIFIED"  >> $GIT_HEADER
echo "#define $PREFIX""_REPO_UNTRACKED  $GIT_UNTRACKED"  >> $GIT_HEADER

if [ $GIT_MODIFIED = 1 ]; then
echo "" >> $GIT_HEADER
echo "#if $PREFIX""_REPO_MODIFIED == 1" >> $GIT_HEADER
echo "  #ifndef "$PREFIX"_GITHEAD_NOWARNING" >> $GIT_HEADER
echo "    #define "$PREFIX"_GITHEAD_NOWARNING" >> $GIT_HEADER
echo "    #pragma message (\"\n\"\\" >> $GIT_HEADER
echo "                     \"******************************************\n\"\\" >> $GIT_HEADER
echo "                     \"*                                        *\n\"\\" >> $GIT_HEADER
centerName "WARNING: repository modified" >> $GIT_HEADER
centerName "changes not committed in"  >> $GIT_HEADER
centerName $REPO_NAME >> $GIT_HEADER
echo "                     \"*                                        *\n\"\\" >> $GIT_HEADER
echo "                     \"******************************************\n\"\\" >> $GIT_HEADER
echo "                    )" >> $GIT_HEADER
echo "  #endif" >> $GIT_HEADER
echo "#endif" >> $GIT_HEADER
fi
if [ $GIT_UNTRACKED = 1 ]; then
echo "" >> $GIT_HEADER
echo "#if $PREFIX""_REPO_UNTRACKED == 1" >> $GIT_HEADER
echo "  #ifndef "$PREFIX"_GITHEAD_NO_UNTRACKED_WARNING" >> $GIT_HEADER
echo "    #define "$PREFIX"_GITHEAD_NO_UNTRACKED_WARNING" >> $GIT_HEADER
echo "    #pragma message (\"\n\"\\" >> $GIT_HEADER
echo "                     \"******************************************\n\"\\" >> $GIT_HEADER
echo "                     \"*                                        *\n\"\\" >> $GIT_HEADER
centerName "WARNING: repository untracked" >> $GIT_HEADER
centerName "untracked files exist in"  >> $GIT_HEADER
centerName $REPO_NAME >> $GIT_HEADER
echo "                     \"*                                        *\n\"\\" >> $GIT_HEADER
echo "                     \"******************************************\n\"\\" >> $GIT_HEADER
echo "                    )" >> $GIT_HEADER
echo "  #endif" >> $GIT_HEADER
echo "#endif" >> $GIT_HEADER
fi

echo "" >> $GIT_HEADER
echo "#endif" >> $GIT_HEADER
