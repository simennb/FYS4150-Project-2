import numpy as np
import matplotlib.pylab as plt

N = [10, 50, 100, 200, 300, 400, 500]
iterations = [70, 3680, 15521, 64025, 145359, 260903, 410155]

p = np.polyfit(N, iterations, 2)

x = np.linspace(0, 500, 1000)
y = p[0]*x**2 + p[1]*x + p[2]

plt.plot(N, iterations, 'ob', label = 'Iteration points')
plt.plot(x, y, '-g', label = r"$1.67x^2 -17.73x + 325.29$")
plt.grid('on')
plt.legend(loc=2)
plt.xlabel('n-values', size = 18)
plt.ylabel('Iterations', size = 18)
plt.title('Iterations as a function of n', size = 16)
plt.savefig('../figures/iterations.png')
plt.show()
