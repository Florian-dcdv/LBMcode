#!/bin/bash

for rad in "20" "30" "40" "50"
do
    if [ -d rad$rad ]; then
        rm -rf rad$rad
    fi
    cp -r defaultCase rad$rad
    cd rad$rad
    sed -i "s=CHARRAD=$rad=g" input.command 
    sed -i "s=CHARRAD=$rad=g" runCase.cmd
    qsub runCase.cmd
    cd ..
done

