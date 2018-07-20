#!/bin/bash
#PBS -e /home/arthbous/atm2d-v0.98b/run_physic4a/test.error
#PBS -o /home/arthbous/atm2d-v0.98b/run_physic4a/test.output
#PBS -d /home/arthbous/atm2d-v0.98b/run_physic4a/

python2.7 ../main.py /home/arthbous/atm2d-v0.98b/run_physic4a
