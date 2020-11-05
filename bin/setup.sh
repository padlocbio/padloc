#!/bin/bash

# script/setup: Set up application for the first time after cloning, or set it
#               back to the initial first unused state.

set -e

cd "$(dirname "$0")/.."

info() {
  printf "$(date "+(%H:%M:%S)") >> ${@}"
}

./bin/bootstrap.sh

./bin/updatedb.sh

printf "\nSet up finished, add padloc to your \$PATH with: \n\n"
printf "    export PATH=\"%s:\$PATH\"\n\n" $(pwd)
printf "Then try typing:\n\n"
printf "    padloc  --help\n"
