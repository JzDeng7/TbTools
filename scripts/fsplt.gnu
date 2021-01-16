set term pdfcairo  lw 1 font "monospace" 
set size square
#"Times New Roman,24" size 4,3 

set yrange [*:*] noextend
set xrange [*:*] noextend
set key off

  set palette defined ( 0 "white", 0.000001 "red", 1 "red")
# set palette defined (-1 "red", -0.81 "red", -0.70 "white", -0.69 "red", 1 "red")
# set palette defined (-1 "red", 0 "white", 1 "red")
# Output
set output "fs.pdf"
  plot "fs.dat" u 1:2:3 with image
unset output
