set size 0.65, 1 
set xlabel "i" # x-axis set ylabel "j" 
set zrange [-400:400]
set cntrparam levels 10 
#set title "E_z" 
#splot "EData5.txt" using 1:2:3 w l 

set term gif animate delay 1
set output "e_movie.gif"

do for [i = 1:19] {
    input = sprintf("EData%d.txt",i)
    splot input using 1:2:3 w l 
}
