#!/bin/bash

# delete output folder 
rm -rf impact_output_moving/ 
# create new output folder
mkdir impact_output_moving   
# everything from -L$BASILISK/gl onward is just for visualisation (but needs a correct installation)
qcc -autolink -O2 -Wall -fopenmp -g -Wall -pipe -D_FORTIFY_SOURCE=2 -o impact_output_moving/3D_Impact 3D_Impact.c -L$BASILISK/gl -lglutils -lfb_tiny -lGLU -lGLEW -lGL -lX11 -lm

# enter folder
cd impact_output_moving  
mkdir Slices
mkdir Interfaces
mkdir Animations

# execute code given resources
export OMP_NUM_THREADS=8

# Parameters are:
# 1: maxlevel 
# 2: impingement angle
# 3: initial drop (vertical) velocity - dimensional, m/s
# 4: pool (horizontal) velocity - dimensional, m/s
# 5: drop radius - dimensional, m
# 6: pool depth - dimensional, m
# 7: domain size non-dimensional, dimensionless, relative to radius
# 8: simulation end time, dimensionless
./3D_Impact 9 90.0 2.45 0.15 1.1e-3 1.65e-3 4.0 0.801

# exit folder      
cd ..           
