#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# Generate copies of a source signal file for multiple seismic shots.
#
# When an external source signal is used, DENISE requires a separate source
# signal file for each shot, even if all shots use the same source wavelet.
# This script automatically generates the required files by copying the input
# source signal.
#
# Usage:
#   ./generate_shot_signals.sh INPUT_SIGNAL NSHOTS OUTPUT_PREFIX
#
# Arguments:
#   INPUT_SIGNAL   Input ASCII source signal (e.g., wavelet file).
#   NSHOTS         Number of shot files to generate.
#   OUTPUT_PREFIX  Prefix for the generated files.
#
# Example:
#   ./generate_shot_signals.sh wavelet.dat 10 source_name
#
# Generates:
#   source_name_shot_1.dat
#   source_name_shot_2.dat
#   ...
#   source_name_shot_10.dat
# -----------------------------------------------------------------------------

if [ $# -ne 3 ]; then
    echo "Usage: $0 INPUT_SIGNAL NSHOTS OUTPUT_PREFIX"
    exit 1
fi

insignal="$1"
nshots="$2"
outsignal="$3"

# Check whether the input source signal file exists.
if [ ! -f "$insignal" ]; then
    echo "Error: input source signal file '$insignal' not found."
    exit 1
fi

# Create one source signal file for each shot.
for (( ishot=1; ishot<=nshots; ishot++ )); do
    new_name="${outsignal}_shot_${ishot}.dat"
    echo "Creating: $new_name"
    cp "$insignal" "$new_name"
done

echo "Done. Generated $nshots source signal files from '$insignal'."