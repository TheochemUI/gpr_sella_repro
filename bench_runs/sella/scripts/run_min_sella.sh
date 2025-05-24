#!/usr/bin/env bash

set -euo pipefail

IS_HPC="${IS_HPC:-false}"

CWD=$(pwd)
FULL_RUNDIR="$CWD/runs/${SCF_THRESH}/sella/${RUNDIR}/${SPIN}/${INDEX}/"

case "$IS_HPC" in
    true|True|TRUE|1|yes|Yes|YES)
        echo "Running on HPC"

        # Get the job ID (example using SLURM)
        JOB_ID="${SLURM_JOB_ID:-$$}"

        # Define the full path to the scratch directory
        SCRATCH_BASE="/scratch/users/rog32"  # **Change this to your actual scratch base path**
        SCRATCH_DIR="$SCRATCH_BASE/eon_scratch_$JOB_ID"
        mkdir -p "$SCRATCH_DIR/data_eon_sella"
        rsync -avv "${DATZIP}" "$SCRATCH_DIR/data_eon_sella/"
        ls -al "$SCRATCH_DIR/data_eon_sella"
        SCRATCHZIP="$SCRATCH_DIR/data_eon_sella/$(basename ${DATZIP})"


        echo "Using scratch directory: $SCRATCH_DIR"
        echo "Using cwd: $CWD"
        echo "Using FULL_RUNDIR: $FULL_RUNDIR"
        echo "Using SCRATCHZIP: $SCRATCHZIP"
        echo "Using DATAZIP: $DATZIP"

        # Create the scratch directory
        mkdir -p "$SCRATCH_DIR/runner"

        # Run in the scratch directory using full paths
        pushd "$SCRATCH_DIR/runner"
        ln -sf ~/Git/Github/TheochemUI/gprd_sella_bench/{.pixi,pixi.toml,pixi.lock} "$SCRATCH_DIR"
        pixi run python "${CWD}/run_sella.py" "$INDEX" --spin "$SPIN" --scfthresh "${SCF_THRESH}" --min_dir "${MIN_DIR}"
        popd

        rsync -a "$SCRATCH_DIR/runner/" "$FULL_RUNDIR/"

        # Clean up the scratch directory
        rm -rf "$SCRATCH_DIR"
        ;;
    false|False|FALSE|0|no|No|NO|"")
        echo "Running locally (IS_HPC=false, 0, no, or empty)"
        cd "${FULL_RUNDIR}"
        pixi run python "${CWD}/run_sella.py" "$INDEX" --spin "$SPIN" --scfthresh "${SCF_THRESH}" --min_dir "${MIN_DIR}"
        ;;
    *)
        echo "Not Running (IS_HPC set to an unexpected value: $IS_HPC)"
        ;;
esac
