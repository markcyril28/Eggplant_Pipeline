#!/bin/bash
set -euo pipefail

FastTree -nt -gtr -gamma -boot 1000 JRO.fas > JRO_boot_1000.newick
