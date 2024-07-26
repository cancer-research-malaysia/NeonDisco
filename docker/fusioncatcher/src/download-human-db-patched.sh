#!/usr/bin/env bash

# download references for Fusioncatcher. Only run once.

if [ -z "$1" ]; then
  echo "No argument provided. Please provide an argument."
  exit 1
else
  $1="/refs"
  # FusionCatcher v1.33
  #
  cd "$(dirname "$(realpath -s "$0")")"
  set -e
  wget --no-verbose --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.aa -O human_v102.tar.gz.aa
  wget --no-verbose --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ab -O human_v102.tar.gz.ab
  wget --no-verbose --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ac -O human_v102.tar.gz.ac
  wget --no-verbose --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.tar.gz.ad -O human_v102.tar.gz.ad
  wget --no-verbose --no-check-certificate http://sourceforge.net/projects/fusioncatcher/files/data/human_v102.md5 -O human_v102.md5
  md5sum -c human_v102.md5

  if [ "$?" -ne "0" ]; then
    echo -e "\n\n\n\033[33;7m   ERROR: The downloaded files from above have errors! MD5 checksums do not match! Please, download them again or re-run this script again!   \033[0m\n"
    exit 1
  fi

  cat human_v102.tar.gz.* > human_v102.tar.gz
  rm -f human_v102.tar.gz.*
  mkdir -p human_v102

  if ! tar -xzf human_v102.tar.gz -C .; then
    echo -e "\n\n\n\033[33;7m   ERROR: The downloaded files are corrupted! Please, download them again or re-run this script again!   \033[0m\n"
    exit 1
  fi

  rm -f human_v102.tar.gz
  rm -f human_v102.md5

  if mv human_v102 /references; then
    "References have been set up. Ready to run Fusioncatcher!"
  else
    "References failed to be set up. Redownload again."
  fi

fi