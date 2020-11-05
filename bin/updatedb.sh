#!/bin/bash

# utils/updatedb: Download and compile the latest version of PADLOC-DB

set -e

cd "$(dirname "$0")/.."

info() {
  printf "$(date "+(%H:%M:%S)") >> ${@}"
}

info "Installing database\n"

[[ -d data ]] && rm -dr data

git submodule update --init --recursive > /dev/null 2>&1

info "Compiling database\n"

find data/hmm/ -name "*.hmm" -exec cat {} > padlocdb.hmm \;
find data/hmm/ -name "*.hmm" -delete
mv padlocdb.hmm data/hmm/
