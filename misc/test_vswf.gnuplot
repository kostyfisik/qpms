f = 'out'
fc = 'outcs'
do for [t=3:137] {
y = ((t-3) % 45)/3
typ = (t-3) % 3
n = floor(sqrt(y+1))
m = y - (n*(n+1)-1)
print 'n = ', n, ', m = ', m, ', typ ', typ
plot f using 1:t w linespoints, fc using 1:t w linespoints
pause -1
}
