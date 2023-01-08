#!/bin/bash 
 
# script file for BASH 
# which bash
# save this file as m.sh
# chmod +x .sh
# ./m.sh
# checked in https://www.shellcheck.net/


printf "compile c file \n"
gcc mu-atom.c -std=c99 -Wall -Wextra -pedantic -O3  -lm

if [ $? -ne 0 ]
then
    echo ERROR: compilation failed !!!!!!
    exit 1
fi

printf "run the compiled program\n"
for k in {1..10}
do
  ./a.out "$k" > mu-atom-"$k".ppm 
done







printf "convert all ppm files to png using Image Magic v 6 convert \n"
# for all ppm files in this directory
for file in *.ppm ; do
  # b is name of file without extension
  b=$(basename "$file" .ppm)
  # convert  using ImageMagic : -resize widthxheight || 
  convert "${b}".ppm "${b}".png  # iWidth = iHeight* DisplayAspectRatio 
  echo "$file"
done


printf "delete all ppm files \n"
rm ./*.ppm


