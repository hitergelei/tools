#!/usr/bin/env python
import numpy as np
import tensorflow as tf

def f(x, a, mu, sigma):
    return a/sigma/np.sqrt(2.*np.pi) *np.exp(-0.5*np.square((x-mu)/sigma))
## User variables
max_step       = 5000
optimizer      = tf.contrib.opt.NadamOptimizer
learnrate      = 1e-3
print_interval = max_step / 10
ngauss         = 2
npz_name       = 'wrapped-gete-conv-md.traj_slice-0-None-1_sym-Ge-Te_nBin-500_rcut-7.1_.npy'

## initial guesses
gmeans  = np.random.rand(ngauss)+2.5
gsigma  = np.random.rand(ngauss)/5.+0.1
gcoeffi = np.random.rand(ngauss)/5.+1.

## Main
curve = np.load(npz_name)
curve = curve[150:300]
x = np.array(curve[:,0])
y_hat = curve[:,1]

# tf variables
mu = tf.Variable(gmeans)
sigma = tf.Variable(gsigma)
a = tf.Variable(gcoeffi)
# function
y = [a[i]/sigma[i]/np.sqrt(2.*np.pi) *tf.exp(-0.5*tf.square((x-mu[int(i)])/sigma[i])) for i in range(ngauss)]
# cost function
cost = tf.norm(tf.reduce_sum(y,axis=0)-y_hat)
# optimizer
optimize = optimizer(learnrate).minimize(cost)

#initialize
init = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init)

print('___step_|______________a______________|______________mu_____________|____________sigma____________|__cost______')
## Train
for i in range(max_step):
    sess.run(optimize)
    if i % print_interval == 0:
        print('{:7d} | {:7.5f}   {:7.5f}           | {:7.5f}   {:7.5f}           | {:7.5f}   {:7.5f}           | {:7.5f}'.format(
        # print('{:5d} | {:7.5f}   {:7.5f}   {:7.5f} | {:7.5f}   {:7.5f}   {:7.5f} | {:7.5f}   {:7.5f}   {:7.5f} | {:7.5f}'.format(
            i, *np.squeeze(sess.run(a)), *np.squeeze(sess.run(mu)), *np.squeeze(sess.run(sigma)), sess.run(cost)))
print('  Final | {:7.5f}   {:7.5f}           | {:7.5f}   {:7.5f}           | {:7.5f}   {:7.5f}           | {:7.5f}'.format(
# print('Final | {:7.5f}   {:7.5f}   {:7.5f} | {:7.5f}   {:7.5f}   {:7.5f} | {:7.5f}   {:7.5f}   {:7.5f} | {:7.5f}'.format(
    *np.squeeze(sess.run(a)), *np.squeeze(sess.run(mu)), *np.squeeze(sess.run(sigma)), sess.run(cost)))

from matplotlib import pyplot as plt
font = {'family':'Arial'}
plt.rc('font', **font)
plt.plot(curve[:,0], curve[:,1], c='k', label='DFT-MD 600K')
plt.plot(curve[:,0], np.sum(sess.run(y),axis=0), c='r', linestyle='dashed', label='$y_{1}+y_{2}$')
plt.plot(curve[:,0], f(curve[:,0], sess.run(a)[0], sess.run(mu)[0], sess.run(sigma)[0]), linestyle='dotted',
    label='$y_{1}$\n'+'a={:.2f}\n$\mu$={:.2f}\n$\sigma$={:.2f}'.format(sess.run(a)[0], sess.run(mu)[0], sess.run(sigma)[0]))
plt.plot(curve[:,0], f(curve[:,0], sess.run(a)[1], sess.run(mu)[1], sess.run(sigma)[1]), linestyle='dotted',
    label='$y_{2}$\n'+'a={:.2f}\n$\mu$={:.2f}\n$\sigma$={:.2f}'.format(sess.run(a)[1], sess.run(mu)[1], sess.run(sigma)[1]))
# plt.plot(curve[:,0], f(curve[:,0], sess.run(a)[2], sess.run(mu)[2], sess.run(sigma)[2]), linestyle='dotted')
plt.grid(alpha=0.3)
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.legend(fontsize='large')
plt.show()
