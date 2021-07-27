#!/bin/bash

# script/bootstrap: Download dependencies required by padloc.

set -e

cd "$(dirname "$0")/.."

info() {
  printf "$(date "+(%H:%M:%S)") >> ${@}"
}

check_command() {
    command -v "$1" > /dev/null 2>&1;
}

info "Installing dependencies\n"

if ! check_command brew; then
	info "Installing Homebrew\n"
	/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
	info "Homebrew installed (%s)\n" $(brew -v | cut -d ' ' -f 2 | head -n 1)
else
	info "Homebrew already installed (%s)\n" $(brew -v | cut -d ' ' -f 2 | head -n 1)
fi

info "Updating Homebrew\n"
brew update 2>&1 >/dev/null

info "Installing Homebrew dependencies\n"
if ! check_command R; then 
	info "Installing R\n"
	brew install R
	info "R installed (%s)\n" $(R --version | cut -d ' ' -f 3 | head -n 1) 
else
	info "R already installed (%s)\n" $(R --version | cut -d ' ' -f 3 | head -n 1) 
fi
if ! check_command hmmsearch; then 
	info "Installing HMMER\n"
	brew install hmmer
	info "HMMER installed (%s)\n" $(hmmsearch -h | grep HMMER | cut -d ' ' -f 3)
else
	info "HMMER already installed (%s)\n" $(hmmsearch -h | grep HMMER | cut -d ' ' -f 3)
fi
if ! check_command prodigal; then 
	info "Installing Prodigal\n"
	brew install prodigal
	info "Prodigal installed (%s)\n" $(prodigal -v 2>&1 >/dev/null | grep -o '[0-9].[0-9].[0-9]') 
else
	info "Prodigal already installed (%s)\n" $(prodigal -v 2>&1 >/dev/null | grep -o '[0-9].[0-9].[0-9]')
fi

info "Installing R dependencies\n"

info "Installing tidyverse\n"
R --slave -e 'install.packages("tidyverse", repos = "http://cran.us.r-project.org", quiet = TRUE)'
info "Installing getopt\n"
R --slave -e 'install.packages("getopt", repos = "http://cran.us.r-project.org", quiet = TRUE)'
info "Installing yaml\n"
R --slave -e 'install.packages("yaml", repos = "http://cran.us.r-project.org", quiet = TRUE)'
