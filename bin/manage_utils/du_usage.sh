#!/bin/bash
# -----------------------------------------------------------------------------
# du_usage.sh — Submit a SLURM job that reports disk usage of each
#               subdirectory, largest first.
#
# Usage:
#   bash du_usage.sh [TARGET_DIR]     # or: chmod +x du_usage.sh && ./du_usage.sh
#
#   TARGET_DIR  Directory to scan (default: current directory).
#               The output file disk_usage_<jobid>.txt is written there too.
# -----------------------------------------------------------------------------

set -eu

TARGET_DIR="${1:-$PWD}"

# Single quotes around --wrap are required: $SLURM_CPUS_PER_TASK must be
# expanded inside the job, not at submission time (where it is unset).
sbatch \
    --account=def-jacquese \
    --job-name=du_usage \
    --time=20:00:00 \
    --mem=2G \
    --cpus-per-task=4 \
    --chdir="$TARGET_DIR" \
    --output=disk_usage_%j.txt \
    --wrap='ls -d -- */ | xargs -P "$SLURM_CPUS_PER_TASK" -I{} du -sh {} 2>/dev/null | sort -rh'
