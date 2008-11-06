#!/bin/bash

make dist
scp -P 2001 RNAcode-0.1.tar.gz washietl@localhost:~/programs/
rsh -p 2001 washietl@localhost /homes/washietl/programs/compileRNAcode.sh
rsh -p 2001 washietl@localhost RNAcode -V

