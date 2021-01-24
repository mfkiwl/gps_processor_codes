q(t, x, dx, n) = ( 1e-15 < abs(dx) && 0 < n ) ? q(t, x+dx, -(x + sin(x)- t)/(1+cos(x)), n-1): x/2

set palette defined (0 "blue",17 "#00ffff",33 "white",50 "yellow",66 "red", 100 "#990000",101 "gray")
set style line 112 lt 0 lw 0.5 lc rgb "#006600"
unset key
unset border
unset tics
set cbtics
set cbrange [-4:8]
set pm3d nohidden3d
set view 0,0,1.6
set colorbox user size 0.02,0.3 orig 0.92,0.65
set term pngc enh font "Arial,8" size 480,240
set out "heatmap.png"

splot \
'temperature.dat' using \
  (th=q(pi*sin($1*pi/180),2*asin($1/90.0),3e-14,10),($2/90.0)*cos(th)): (sin(th)):(1):($3)  w pm3d,\
'world.dat' using \
  (th=q(pi*sin($2/180*pi),2*asin($2/90),3e-14,10),($1/90)*cos(th)): (sin(th)):(2)  w lines lc rgb "black"

! display heatmap.png