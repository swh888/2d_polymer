#set encoding iso_8860_1
#set terminal postscript solid color enhanced "Helvetica, 30"
set terminal png enhanced

stats 'conf.txt'
NMAX=STATS_blocks-2

#Set labels
#set xlabel "{/Helvetica-Oblique x}" font ",42" offset 1
#set ylabel "{/Helvetica-Oblique y}" font ",42" offset 0.6
unset xlabel
unset ylabel


#Set tics
set xtics 50
set xtics font "" 
set xtics mirror
set ytics 50
set ytics font ""
set ytics mirror
set size ratio -1
set xrange [-0.5:320]
set yrange [-0.5:320]
#set key font ",22"
unset key
#Plot

do for [n=0:NMAX] {
    outfile=sprintf("rep-%d.png", n)
    set output outfile
    plot 'conf.txt' i n u 2:3 w l lw 2 lc rgb "black"
}

