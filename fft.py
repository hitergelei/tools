#!/usr/bin/env python

import numpy as np
from numpy import fft

## Global params
# Time
t_init = 0
t_fin  = 100
d_t    = 0.1
# f(t)
def ftn(t):
    # return np.exp(-1 * t**2 /2 ) 
    return np.exp(1.j * t)
    # return np.cos(10.*t)
# Pyplot
xlabel_1 = 'Time (s)'
xlabel_2 = 'Frequency (Hz)'

#===================== Main code ====================
"""
 Follow Arfken conventions
            1           -inf                              
y(f) = ------------ * integral {y(t) * exp(2*pi*f*t) * dt}
        sqrt(2*pi)       inf                              

        e.g) 1/sqrt(2*pi) exp(iat) --> delta(f-a)
        ref) https://www.chegg.com/homework-help/questions-and-answers/1-list-fourier-transform-pairs-table-41-p-252-f-t-even-function-observe-nature-fourier-tra-q11614127
"""
## Calculation
# dtype management
d_t = float(d_t)
# Get discrete time domain
t = np.arange(t_init, t_fin + d_t, d_t)
# Get frequency domain (Use Arfken representation)
f = fft.fftfreq(len(t)) / d_t
f = fft.fftshift(f)

# Get f(t)
y_t = ftn(t)
# Get f(f) and normalize
y_f = fft.fft(y_t, norm='ortho') * np.sqrt(d_t * (t_fin - t_init + d_t) / 2 / np.pi)
y_f = fft.fftshift(y_f)

# Test options
# print(np.trapz(((y_f)),f))
# print(f[-1])
# print(y_f[len(t)/2])

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
ax[1,0].plot(f, np.abs(y_f))
ax[1,0].set_xlabel(xlabel_2)
ax[1,0].set_ylabel(ylabel_ampli)
ax[1,0].grid(True)
# Plot transformed function's phase
ang = np.angle(y_f)
ax[1,1].plot(f, ang * (ang >= 0) + (ang+2*np.pi) * (ang < 0))
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
ax[3,0].plot(f, np.real(y_f))
ax[3,0].set_xlabel(xlabel_2)
ax[3,0].set_ylabel(ylabel_real)
ax[3,0].grid(True)
# Plot transformed function's imaginary part
ax[3,1].plot(f, np.imag(y_f))
ax[3,1].set_xlabel(xlabel_2)
ax[3,1].set_ylabel(ylabel_imag)
ax[3,1].grid(True)

# Plot
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=1.0)
plt.show()

