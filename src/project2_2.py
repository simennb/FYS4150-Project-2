from pylab import *
import numpy as np
import matplotlib.pyplot as plt
rc('font',**{'family':'serif'}) 

n = 200

interact = ['two-noint', 'two-int']
colors = ['blue','black']
lines = ['-','--']
plt.figure()
ax = plt.gca()

for k in range(2):

    psi = zeros((4,n))
    omega_r = []

    filename = '../benchmarks/eigenvectors_%s_n%d.dat'%(interact[k],n)
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

        for i in range(len(omega_r)-2):  # take only two first for readability
            plot(rho,abs((psi[i,:])**2),linewidth=1.4, ls='%s'%lines[k],color='%s'%colors[i])
title(r'Wavefunction of gs for different $\omega_r$, $N=$%d'%n,size=16)
xlabel(r'$\rho$',size=18)
ylabel(r'$|\Psi(\rho)|^2$',size=18)

legend([r'$\omega_r =$ %.2f'%omega_r[0],r'$\omega_r =$ %.2f'%omega_r[1]],fontsize=15)
leg = ax.get_legend()
leg.legendHandles[0].set_color('blue')
leg.legendHandles[1].set_color('black')
grid('on')
xlim([0,4])
savefig('../figures/eigvec_compare_n%d.png'%n)
show()
