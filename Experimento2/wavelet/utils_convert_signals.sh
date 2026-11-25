#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# signal_utils.sh
#
# Utilities for converting and visualizing source signals used by DENISE.
#
# DENISE expects source signals in ASCII format, while Seismic Unix (SU)
# typically works with binary float and SU trace formats. This script provides
# a quick reference for the most common conversion operations.
# -----------------------------------------------------------------------------

echo "=============================================="
echo " Signal utilities for DENISE and Seismic Unix "
echo "=============================================="
echo
echo "1. ASCII (.dat) -> Binary float (.bin)"
echo "--------------------------------------"
echo "a2b < signal.dat > signal.bin outpar=/dev/null"
echo
echo "Optional:"
echo "  n1=<samples_per_line>   (default: 2)"
echo
echo "2. Binary float (.bin) -> SU (.su)"
echo "----------------------------------"
echo "suaddhead ns=<NS> < signal.bin | \\"
echo "sushw key=dt,tracf a=<DT>,1 > signal.su"
echo
echo "3. Display SU trace"
echo "-------------------"
echo "suxwigb < signal.su"
echo
echo "or"
echo
echo "suximage < signal.su"
echo
echo "4. Display header information"
echo "-----------------------------"
echo "surange < signal.su"
echo
echo "5. SU (.su) -> Binary float (.bin)"
echo "----------------------------------"
echo "sustrip < signal.su > signal.bin"
echo
echo "6. Binary float (.bin) -> ASCII (.dat)"
echo "--------------------------------------"
echo "b2a < signal.bin > signal.dat"
echo
echo "=============================================="