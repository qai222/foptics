#!/usr/bin/env bash

awk 'BEGIN{i=0} /<dielectricfunction>/,\
                /<\/dielectricfunction>/ \
                 {if ($1=="<r>") {a[i]=$2 ; xx[i]=$3 ; yy[i]=$4 ; zz[i]=$5 ; xy[i]=$6; yz[i]=$7; xz[i]=$8; i=i+1}} \
     END{for (j=0;j<i/2;j++) print a[j], xx[j], xy[j], xz[j], xy[j], yy[j], yz[j], xz[j], yz[j], zz[j]}' vasprun.xml > im.epsilon


awk 'BEGIN{i=0} /<dielectricfunction>/,\
                /<\/dielectricfunction>/ \
                 {if ($1=="<r>") {a[i]=$2 ; xx[i]=$3 ; yy[i]=$4 ; zz[i]=$5 ; xy[i]=$6; yz[i]=$7; xz[i]=$8; i=i+1}} \
     END{for (j=i/2;j<i;j++) print a[j], xx[j], xy[j], xz[j], xy[j], yy[j], yz[j], xz[j], yz[j], zz[j]}' vasprun.xml > re.epsilon
