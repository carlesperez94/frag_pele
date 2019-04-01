#!/bin/bash
#BSUB -J "AdaptiveTest"
#BSUB -o AdaptiveTest.out
#BSUB -e AdaptiveTest.err
#BSUB -n 5
#BSUB -W 00:50
#BSUB -q bsc_debug


python runAllTests.py
