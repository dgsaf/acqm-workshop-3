#
set terminal epslatex input color solid

#
set xlabel "$x$ [$a_{0}$]"
set ylabel "Energy [Ha]"
set xrange [-5:5]
set yrange [ -1:5]
unset key

#
wf(nx,ty)=sprintf("../../output/qho/wf.n_x-%i.%s.txt", nx, ty)
v(x) = 0.5*(x**2)
k = 15

# shooting bisection
set output "wf_sb.tex"
set title "Vibrational Wavefunctions [Shooting-Bisection]"

plot v(x) lt rgb "black", \
  for [n=0:3] wf(2**k,"sb") u 1:(0.5*column(n+2)+n+0.55) \
    w lines lt rgb "red" , \
  for [n=0:3] wf(2**k,"an") u 1:(0.5*column(n+2)+n+0.5) \
    w lines lt rgb "blue"

# numerov-cooley
set output "wf_nc.tex"
set title "Vibrational Wavefunctions [Numerov-Cooley]"

plot v(x) lt rgb "black", \
  for [n=0:3] wf(2**k,"nc") u 1:(0.5*column(n+2)+n+0.55) \
    w lines lt rgb "red" , \
  for [n=0:3] wf(2**k,"an") u 1:(0.5*column(n+2)+n+0.5) \
    w lines lt rgb "blue"
