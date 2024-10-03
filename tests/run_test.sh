#!/usr/bin/env bash

set -e

if [ -d tmpout ]
then
    rm -r tmpout
fi
mkdir tmpout

python3 python/test_SM.py
python3 python/test_SMEFT.py
./bash/test_cfgs.sh
./bash/test_fit.sh
