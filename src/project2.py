from pylab import *
import numpy as np
rc('font',**{'family':'serif'}) 

n = 100

psi = zeros((4,n))
omega_r = []

filename = '../benchmarks/eigenvectors_int_n%d.dat'%n
data = open(filename,'r')
count1 = 0
count2 = 0
count_omega = -1
for line in data:
    if count1 == 0:
        nuline = line.split(' ')
        n = int(nuline[3].strip())
        rho_min = float((nuline[5]))
        rho_max = float((nuline[-1])[:-2])
    if line.startswith('omega_r'):
        omega_r.append(float((line.split())[1]))
        count2 = 0
        count_omega += 1
    elif count1 != 0:
        psi[count_omega, count2] = float(line)
        count2 += 1
    count1 += 1

rho = linspace(rho_min,rho_max,n) # never stored it, but equally spaced anyway

for i in range(len(omega_r)):
    plot(rho,abs((psi[i,:])**2), label=r'$\omega_r =$%.2f'%omega_r[i],linewidth=2.0)
title(r'Wavefunction of gs for different $\omega_r$, $n=$%d'%n,size=16)
xlabel(r'$\rho$',size=18)
ylabel(r'$|\Psi(\rho)|^2$',size=18)
legend(fontsize=15)
grid('on')
savefig('../figures/eigvec_interact_n%d.png'%n)
show()
