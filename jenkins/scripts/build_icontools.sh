#!/bin/bash

set -e -x

function error {
    echo "*** Error: $@" >&2
    exit 1
}

# Check if script is called correctly
[[ $(git rev-parse --show-toplevel 2>/dev/null) = $(pwd) ]] || error "$0 not launched from toplevel of repository"

pushd ext
# Activate spack
. spack-c2sm/setup-env.sh

# Build icontools
spack install icontools@c2sm-master%gcc

