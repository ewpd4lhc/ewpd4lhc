#!/usr/bin/env bash

set -e

../ewpd4lhc.py --output tmpout/testfit.yml --root_output tmpout/testfit.root
../ROOTfit/fit.py -i tmpout/testfit.root -p cHu -s=-0.2:0.2
../ROOTfit/fit.py -i tmpout/testfit.root -p cHu -s=-0.2:0.2 -f cHj1 --plot
../ROOTfit/fit.py -i tmpout/testfit.root -p cHu,cHj1 -s=-0.2:0.2,-0.2:0.2  --plot
