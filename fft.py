#!/usr/bin/env python

import numpy as np
from numpy import fft

## Global params
# Time
t_init = -100
t_fin  = 100
d_t    = 0.1
# f(t)
def ftn(t):
    return np.exp(1.j*t)
# Pyplot
xlabel_1 = 'Time (s)'
xlabel_2 = 'Angular velocity (rad/s)'

#===================== Main code ====================
## Calculation
# dtype management
d_t = float(d_t)
# Get discrete time domain
t = np.arange(t_init, t_fin + d_t, d_t)
# Get angular velocity domain (Use Arfken representation)
w = fft.fftfreq(len(t)) * 2*np.pi / d_t
w = np.concatenate((w[int(round(len(t)/2.)):], w[:int(round(len(t)/2.))]))

# Get f(t)
y_t = ftn(t)
# Get f(w) and normalize
y_w = fft.fft(y_t) / (t_fin - t_init)
y_w = np.concatenate((y_w[int(round(len(t)/2.)):], y_w[:int(round(len(t)/2.))]))

## Plot 
from matplotlib import pyplot as plt
fig, ax = plt.subplots(4,2)
ylabel_ampli = 'Amplitude'
ylabel_phase = 'Phase'
ylabel_real  = 'Re(f)'
ylabel_imag  = 'Im(f)'
# Plot original function's amplitude
ax[0,0].plot(t, np.abs(y_t))
ax[0,0].set_xlabel(xlabel_1)
ax[0,0].set_ylabel(ylabel_ampli)
ax[0,0].grid(True)
# Plot original function's phase
ang = np.angle(y_t)
ax[0,1].plot(t, ang * (ang >= 0) + (ang+2*np.pi) * (ang < 0))
ax[0,1].set_xlabel(xlabel_1)
ax[0,1].set_ylabel(ylabel_phase)
ax[0,1].set_ylim(0,2*np.pi)
ax[0,1].grid(True)
# Plot transformed function's amplitude
ax[1,0].plot(w, np.abs(y_w))
ax[1,0].set_xlabel(xlabel_2)
ax[1,0].set_ylabel(ylabel_ampli)
ax[1,0].grid(True)
# Plot transformed function's phase
ang = np.angle(y_w)
ax[1,1].plot(w, ang * (ang >= 0) + (ang+2*np.pi) * (ang < 0))
ax[1,1].set_xlabel(xlabel_2)
ax[1,1].set_ylabel(ylabel_phase)
ax[1,1].set_ylim(0,2*np.pi)
ax[1,1].grid(True)
# Plot original function's real part
ax[2,0].plot(t, np.real(y_t))
ax[2,0].set_xlabel(xlabel_1)
ax[2,0].set_ylabel(ylabel_real)
ax[2,0].grid(True)
# Plot original function's imaginary part
ax[2,1].plot(t, np.imag(y_t))
ax[2,1].set_xlabel(xlabel_1)
ax[2,1].set_ylabel(ylabel_imag)
ax[2,1].grid(True)
# Plot transformed function's real part
ax[3,0].plot(w, np.real(y_w))
ax[3,0].set_xlabel(xlabel_2)
ax[3,0].set_ylabel(ylabel_real)
ax[3,0].grid(True)
# Plot transformed function's imaginary part
ax[3,1].plot(w, np.imag(y_w))
ax[3,1].set_xlabel(xlabel_2)
ax[3,1].set_ylabel(ylabel_imag)
ax[3,1].grid(True)

# Plot
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=1.0)
plt.show()

