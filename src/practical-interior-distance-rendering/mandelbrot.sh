#!/bin/bash
make
for i in $(seq -w 1 10)
do
  ./mandelbrot 0 $i 1024 1024 -0.75 0 1.5 period-$i.ppm
done
while read tag maxiters view
do
  echo $tag
  echo -n "plain"
  time ./mandelbrot 0 $maxiters 1024 1024 $view ${tag}-0-plain.ppm
   echo -n "unbiased"
  time ./mandelbrot 1 $maxiters 1024 1024 $view ${tag}-1-unbiased.ppm
  echo -n "biased"
  time ./mandelbrot 2 $maxiters 1024 1024 $view ${tag}-2-biased.ppm
  ./mandelbrot 6 $maxiters 1024 1024 $view ${tag}-6-analysis.ppm
  echo
done <<EOF
full 4096 -0.75 0 1.5
ejs 4096 -1.76897871630668612705e+00 4.605239642778105129505e-03 2.1579186437577746e-05
mini 65536 -1.768978682440136083993913575455745429e+00 4.605397058461370118234574128620344171e-03 2.1073424255447021e-08
EOF
