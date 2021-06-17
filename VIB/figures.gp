#
#set terminal epslatex input color solid

# files
ssg=sprintf("../PECS/PEC.1ssg")
psu=sprintf("../PECS/PEC.2psu")

ssg_wf=sprintf("../output/vib/1ssg.txt")
psu_wf=sprintf("../output/vib/2psu.txt")

fc_ssg=sprintf("../output/vib/fc_1ssg.txt")

# parameters
n_wf = 60
n_wf_shown = 10
step = n_wf / n_wf_shown

# common settings
set palette defined (0 "blue" , 1 "red")
unset colorbox
set grid xtics ytics ztics
set key \
  top right \
  box opaque \
  samplen 1 spacing 0.6 height +0.6

# dissociative wavefunction settings
set xlabel "$R$"
set ylabel "Energy [Ha]"

set xrange [0:10]
set yrange [-0.6:1.0]

set format x "\\scriptsize %g"
set format y "\\scriptsize %g"
set key width -3.5 spacing 0.8 height +0.6


# figure(s): 1ssg dissociative
# set output sprintf("figure_1ssg.tex")
unset key
set title sprintf("$1s\\sigma_{g}$ Dissociative Wavefunctions")

plot \
  ssg using 1:2 \
    with lines lc "black" , \
  for [i=2:n_wf:step] ssg_wf using 1:i \
    with lines palette frac ((n_wf-i)/(n_wf-2.0))

pause -1
# set output


# figure(s): 2psu dissociative
# set output sprintf("figure_2psu.tex")
unset key
set title sprintf("$2p\\sigma_{u}$ Dissociative Wavefunctions")

plot \
  psu using 1:2 \
    with lines lc "black" , \
  for [i=2:n_wf:step] psu_wf using 1:i \
    with lines palette frac ((n_wf-i)/(n_wf-2.0))

pause -1
# set output


# franck-condon settings
set xlabel "$Kinetic Energy Release [eV]$"
set ylabel ""

set xrange [0:30]
set autoscale y

set format x "\\scriptsize %g"
set format y "\\scriptsize %g"
set key width -3.5 spacing 0.8 height +0.6

# figure(s): franck-condon
# set output sprintf("figure_fc_1ssg.tex")
unset key
set title sprintf("Franck-Condon Approximation for $1s\\sigma_{g}$")

plot \
  for [i=2:10:3] fc_ssg using 1:i \
    title sprintf("\\nu_{%i}", i-1) \
    with lines palette frac ((10.0-i)/(10.0-2.0))

pause -1
# set output
