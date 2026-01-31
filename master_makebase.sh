#!/bin/bash

source config.sh
source ./code/syntest_toAv.sh
source ./code/syntest_toAv_combinechunks.sh

pycodepath="./code/"
sourcepath="./sourcedata/"

# Generate synthetic data in chunks for both scenarios
syntest_toAv "nomiss" $sourcepath $derivativepath $pycodepath
syntest_toAv "yesmiss" $sourcepath $derivativepath $pycodepath

# Combine chunks into final datasets
echo "Combining chunks for nomiss scenario..."
syntest_toAv_combinechunks "nomiss"

echo "Combining chunks for yesmiss scenario..."
syntest_toAv_combinechunks "yesmiss"

