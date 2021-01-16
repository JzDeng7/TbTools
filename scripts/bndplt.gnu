set term pdfcairo  lw 2 font "monospace" size 4,3
#"Times New Roman,24" size 4,3 

set xrange [0:*] noextend
set key off
# Fermi energy
set arrow from graph 0, first 0 to graph 1, first 0  nohead dt 2 lc "grey"
# Output
set output "bnd.pdf"
  plot "bnd.dat" u 1:2:(column(-2)%7+1) w l lc var
unset output
