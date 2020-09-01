#!/bin/sh
#tput setaf 3; echo "---------------------- COMPILATION ----------------------"; tput sgr0;
make
#tput setaf 2; echo "------------------- PROGRAM EXECUTION -------------------"; tput sgr0;
#echo ""
time ./e
#echo ""
#tput setaf 3; echo "----------------------- CLEANING ------------------------"; tput sgr0;
#make clean