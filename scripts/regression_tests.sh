#! /bin/bash

echo "Starting regression tests"
cd ../tests
ls -l
../src/sMQM -T=0 -V 

echo "Finalized regression tests"
exit 0
