#!/bin/bash

SEQ1_NAME=$1
SEQ2_NAME=$2
XLEN=$3
YLEN=$4
RATIO=$5
FDATA=$6
FIG=$7

gnuplot -e "
 set x2range [0 : ${XLEN}];
 set yrange [${YLEN} : 0];
 set x2tics ${XLEN}<10000 ? 1000 : (${XLEN}<100000 ? 5000 : 10000);
 set ytics  ${YLEN}<10000 ? 1000 : (${YLEN}<100000 ? 5000 : 10000);
 set mx2tics 5;
 set mytics 5;
 set tics out;
 set tics scale 1.5,0.5;
 set x2tics offset 0,-0.5;
 set ytics offset 0.8,0;
 set tics font 'DejaVuLGCSans,7';
 set x2label '${SEQ1_NAME}';
 set x2label font 'DejaVuLGCSans,10';
 set x2label offset 0,-0.8;
 set ylabel '${SEQ2_NAME}';
 set ylabel font 'DejaVuLGCSans,10';
 set ylabel offset 3.3,0;
 set nokey;
 unset xtics;
 set size ratio ${RATIO};
 set terminal png;
 set output '${FIG}';
 plot '${FDATA}' w p ps 0.15 lc rgb 'black'
 "
