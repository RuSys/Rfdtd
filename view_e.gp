set size 0.65, 1 
set xlabel "i" # x-axis set ylabel "j" 
set zrange [-400:400]
set cntrparam levels 10 
#set title "E_z" 
splot "EData8.txt" using 1:2:3 w l 