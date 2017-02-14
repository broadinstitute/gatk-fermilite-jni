#!/usr/bin/env bash

#Helper script to be called from build_both_dylib_and_dylib.sh
#
#checkout and clone a new copy of gatk in a tmpdir
#build the .so file and return its location
#this can only be run on broad machines
#
# usage build_pair_hmm_so_in_clean_repo.sh <commit hash>

set -e
set -v

COMMIT=$1
PROJECT=$2

LIB_PATH="build/libs/jni"

export TMPDIR="/broad/hptmp"

PROJECT_TMP_DIR=`mktemp --tmpdir -d 2>/dev/null || mktemp -d -t 'mytmpdir'`
function finish {
  set +e 
#  rm -rf "$PROJECT_TMP_DIR/$PROJECT"
}
trap finish EXIT

cd "$PROJECT_TMP_DIR"
GIT_LFS_SKIP_SMUDGE=1 git clone git@github.com:broadinstitute/${PROJECT}.git 1>2
cd $PROJECT
GIT_LFS_SKIP_SMUDGE=1 git checkout -f "$COMMIT" 1>2

./gradlew build 1>2

LIB_NAME=$( ls $LIB_PATH ) 
cp "$LIB_PATH/$LIB_NAME" "$PROJECT_TMP_DIR/$LIB_NAME"

echo "$PROJECT_TMP_DIR/$LIB_NAME"
exit 0
