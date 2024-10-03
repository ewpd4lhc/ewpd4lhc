#!/usr/bin/env bash

set -e

for cfg in config/*cfg
do
    echo ">>>>>> Now running $cfg"
    ../ewpd4lhc.py --config $cfg --output tmpout/$(basename ${cfg}).yml
done
