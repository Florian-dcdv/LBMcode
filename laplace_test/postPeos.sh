#!/bin/bash

echo "# rad \t Pliq \t Pgas" > laplaceTest.dat

for rad in "20" "30" "40" "50"
do
    cd rad$rad
    echo "Case: rad$rad"
    Pliq=$(sed -n '101p' ./Lines/t35000_x=0-199y=100-100.dat | awk '{print $7}')
    Pgas=$(sed -n '2p' ./Lines/t35000_x=0-199y=100-100.dat | awk '{print $7}')
    echo $rad $Pliq $Pgas >> ../laplaceTest.dat
    cd ..
done

cat laplaceTest.dat

