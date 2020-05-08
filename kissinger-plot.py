#!/usr/bin/env python

import numpy as np

pairs = [
    # [peak temperature (K), Heat rate (K/s)]
    [780, 7.5e11   ],
    [740, 3.75e11  ],
    [698, 1.875e11 ],
    ]

# > Main
# Data
pairs = np.array(pairs)
x = 1000/pairs[:,0]
y = np.log(pairs[:,1] / pairs[:,0]**2)

# Least squares fit
from sklearn.linear_model import LinearRegression
lr = LinearRegression()
lr.fit(x.reshape(-1,1), y)
slope = lr.coef_[0]
inter = lr.intercept_
from ase.units import kB
Ea = -kB *1000 *slope 

# > Plot
from matplotlib import pyplot as plt
plt.scatter(x, y, c='k')
x_range = np.array([np.min(x), np.max(x)])
plt.plot(x_range, slope*x_range+inter, c='r', label='$E_a$={:.4f}eV'.format(Ea))
plt.title('Kissinger Plot', fontsize='x-large')
plt.xlabel('1000/$T_{p}$ ($K^{-1}$)', fontsize='x-large')
plt.ylabel('ln($\phi/T^{2}_{p}$) ($K^{-1}s^{-1}$)', fontsize='x-large')
plt.legend(fontsize='x-large')
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.subplots_adjust(left=0.16, bottom=0.15, right=0.94, top=0.91)
plt.grid(alpha=0.2)
plt.show()
