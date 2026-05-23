#! /bin/bash
set -euo pipefail

DEB="mega-cc_12.1.2-1_amd64.deb"

sudo apt-get update
sudo apt-get install -y ./"$DEB"
