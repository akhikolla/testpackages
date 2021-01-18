#!/bin/bash
export OMP_NUM_THREADS=1
R -f Binary/binary.R
rc=$?
if [ $rc != 0 ]; then
printf "${RED}Binary test crashed, abort...${NC}\n"
exit 1
fi

R -f Categorical/categorical.R
rc=$?
if [ $rc != 0 ]; then
printf "${RED}Categorical test crashed, abort...${NC}\n"
exit 1
fi

R -f Contingency/contingency.R
rc=$?
if [ $rc != 0 ]; then
printf "${RED}Contingency test crashed, abort...${NC}\n"
exit 1
fi

R -f Continuous/continuous.R
rc=$?
if [ $rc != 0 ]; then
printf "${RED}Continuous test crashed, abort...${NC}\n"
exit 1
fi

