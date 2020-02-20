#!/bin/sh

for i in {00..29}
do
    rm -rf job_$i
    mkdir job_$i
    ase convert -n $i::30 cubic-not-vacuum.traj job_$i/inits.traj
    cp run.sh INCAR POTCAR sub_routine.py job_$i/
    sed "s/ran_gst_liquid/job_$i/" "run.sh" > job_$i/run.sh
    cd job_$i
    qsub run.sh
    cd ..
done
