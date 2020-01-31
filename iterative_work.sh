#!/bin/sh

for i in {02..20}
do
    rm -rf job_$i
    mkdir job_$i
    cp 1-training-set-making.py 1-training-set-making.sh INCAR POTCAR cubic-backbone-64.vasp job_$i
    sed "s/ran_gst_liquid/job_$i/" "1-training-set-making.sh" > job_$i/1-training-set-making.sh
    qsub job_$i/1-training-set-making.sh
done
