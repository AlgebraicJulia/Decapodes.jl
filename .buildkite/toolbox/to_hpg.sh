#!/bin/bash

if [ "$1" == "-t" ]; then
    echo "Flag -t is provided."
    DRYRUN="--dry-run"
  else
    DRYRUN=""
fi

# check if we are in a valid git repository 
REPO=`git rev-parse --show-toplevel`
if [[ ! $REPO ]]; then 
  echo "We are not in a valid git repo $(pwd)"  
  exit 1
fi

REMOTE_DIR="$REPO/.buildkite/toolbox/"
CONFIG="$REPO/.buildkite/config.json"

# obtain username
if [ ! -f "$REMOTE_DIR/username" ]; then
  USERNAME=$USER
  read -p "Is \"${USERNAME}\" your HPG username? [y/N]" -n 1 -r -s
  echo

  if [[ ! $REPLY =~ ^[Yy]$ ]]
  then
    read -p "Please provide the username to your HPG account: " -r USERNAME
  fi
else
  USERNAME=`head -1 $REMOTE_DIR/username`
fi

# get branch name
BRANCH=`echo $REPO | xargs basename`

TIME=`date +%s`
rsync -a --exclude '.buildkite/remote/to_hpg.sh' --exclude '.git' --exclude-from '.gitignore' --rsync-path="mkdir -p builds/ && rsync -arv " $REPO $USERNAME@hpg.rc.ufl.edu:~/builds/



