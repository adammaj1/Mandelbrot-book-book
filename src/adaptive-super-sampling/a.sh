#!/bin/bash 
 
# script file for BASH 
# which bash
# save this file as a.sh
# chmod +x a.sh
# ./a.sh
# should be checked in https://www.shellcheck.net/
# 
# https://mathr.co.uk/blog/2014-11-22_adaptive_supersampling_using_distance_estimate.html


printf "compile main program \n"
gcc a.c -std=c99 -Wall -Wextra -pedantic -O3 -fopenmp -lm


printf "display OMP info \n"
export  OMP_DISPLAY_ENV="TRUE"
export  OMP_DISPLAY_ENV="FALSE"

printf "run the compiled program and make pgm files \n"
for depth in 0 1 2 3 4
do
  time ./a.out 960 540 \
  -1.769384146357043e+00 4.230589061544685e-03 1.1169192421368613e-05 \
  1000000 ${depth} > sparse-${depth}.pgm
done
for depth in 0 1 2 3 4
do
  time ./a.out 960 540 \
  -1.747370489 4.733564e-3 3e-7 \
  1000000 ${depth} > dense-${depth}.pgm
done



printf "change Image Magic settings\n"
export MAGICK_WIDTH_LIMIT=100MP
export MAGICK_HEIGHT_LIMIT=100MP

printf "convert all pgm files to png using Image Magic v 6 convert \n"
# for all pgm files in this directory
for file in *.pgm ; do
  # b is name of file without extension
  b=$(basename "$file" .pgm)
  # convert  using ImageMagic
  # convert "${b}".pgm -resize 2000x2000 "${b}".png
   convert "${b}".pgm "${b}".png
  echo "$file"
done


printf "delete all pgm files \n"
rm ./*.pgm

 
echo OK

printf "info about software \n"
bash --version
make -v
gcc --version
convert -version
convert -list resource
# end
