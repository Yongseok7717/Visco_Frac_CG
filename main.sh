#!/bin/bash

###################################################################################
#   If you want to run this main.sh, type and enter                               #
#                                                                                 #
#   chmod u+x ./main.sh                                                           #
#   ./main.sh                                                                     #                    
#                                                                                 #
#   on your terminal.                                                             #
#                                                                                 #
###################################################################################


chmod u+x ./homo_smooth_cg_modi.py
chmod u+x ./homo_nonsmooth_cg.py

#Nonsmooth case with linear polynomial basis
date
./homo_nonsmooth_cg.py -k 1 -i 1 -I 8 -j 3 -J 10

#Nonsmooth case with quadratic polynomial basis
date
./homo_nonsmooth_cg.py -k 2 -i 1 -I 8    -j 3 -J 10

#Smooth case with linear polynomial basis
date
./homo_smooth_cg_modi.py -k 1 -i 1 -I 8 -j 3 -J 10

#Smooth case with quadratic polynomial basis
date
./homo_smooth_cg_modi.py -k 2 -i 1 -I 8 -j 3 -J 10