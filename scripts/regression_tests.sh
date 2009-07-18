#! /bin/bash

echo "Starting regression tests"
cd ../test
ls -l
../src/sMQM -T=0 -V > ../test/MQM_test0.txt

echo "Finalized regression tests"
exit 0
