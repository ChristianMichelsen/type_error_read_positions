#!/usr/bin/env bash
#OUTPREFIX="./out/"
#mkdir -p $(dirname $OUTPREFIX)

#gargammel.pl -n 1000 -se --comp 0,0,1 -damage 0.03,0.4,0.01,0.3 -o ${OUTPREFIX} data/
gargammel.pl -n 1000000 -se --comp 0,0,1 -damage 0.03,0.4,0.01,0.3 -o ./out/ data/
## -n 1000 is the number of reads simulated
## -se -> singl-end. for paired-end just remove this argument
## --comp 0,0,1 sample p(exogenous|microbial) == 0, p(contamination) == 0, p(hostDNA) == 1
## -damage 0.03,0.4,0.01,0.3 v: nick frequency l: length of overhanging ends (geometric parameter) d: prob. of deamination of Cs in double-stranded parts s: prob. of deamination of Cs in single-stranded parts
