#!/bin/bash

nmax=`awk 'BEGIN{a=0}{if(NF==0) a+=1}END{print a/2}' conf.txt`
gnuplot plot-polymer.gnu

#for ((i=0; i<$nmax; i++)); do
    #psconvert rep-${i}
#    convert -density 150 -quality 80 rep-${i}.pdf rep-${i}.jpg
#done
convert -delay 20 -loop 0 *.png -resize 90% rep-animate.gif

#rm rep*.pdf rep*.png rep*.ps
rm rep*.png
