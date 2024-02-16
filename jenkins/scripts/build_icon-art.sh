#!/bin/bash

set -e -x

function error {
    echo "*** Error: $@" >&2
    exit 1
}

# Check if script is called correctly
[[ $(git rev-parse --show-toplevel 2>/dev/null) = $(pwd) ]] || error "$0 not launched from toplevel of repository"

BRANCH=art
GIT_REMOTE=git@github.com:C2SM/icon.git

pushd ext

# Remove icon folder (if existing)
rm -fr icon-art

# Clone icon
git clone --depth 1 --recurse-submodules -b ${BRANCH} ${GIT_REMOTE} icon-art

pushd icon

if [[ $(hostname) == eu-* ]]; then
    ./jenkins/scripts/jenkins_euler.sh -b -fc gcc --configure euler.cpu.gcc.O2
else
    SPACK_TAG=`cat icon/config/cscs/SPACK_TAG`
    . spack-c2sm/setup-env.sh
    spack env activate -p -d config/cscs/spack/${SPACK_TAG}/daint_cpu_nvhpc
    spack install -u build
fi

popd

popd
